#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import generic_filter
import astropy.io.fits as fits
from scipy.stats import linregress

# Function to calculate local standard deviation using a generic filter
def local_std(image, neighborhood_radius):
    def std_filter(values):
        return np.std(values)
    return generic_filter(image, std_filter, size=(2*neighborhood_radius+1, 2*neighborhood_radius+1))

# Function to estimate background noise in a region
def estimate_background_noise(image_data, x_start, y_start, width=10, height=10):
    x_end = x_start + width
    y_end = y_start + height
    background_pixels = image_data[y_start:y_end, x_start:x_end]
    background_noise = np.mean(background_pixels)
    return background_noise

# Function to calculate signal-to-noise ratio (SNR)
def calculate_snr(image_data, background_noise=0, neighborhood_radius=2):
    signal = image_data - background_noise
    noise = local_std(image_data, neighborhood_radius)
    snr = np.divide(signal, noise, out=np.zeros_like(signal), where=noise != 0)
    return snr

# Directory and file names
file_name_888 = 'concolved_887.491MHz_18arcsec_30dorB.fits'
file_name_1367 = 'concolved_1367.49MHz_18arcsec_30dorB.fits'
file_name_1419 = 'concolved_1419.5MHz_18arcsec_30dorB.fits'

# Read image data from FITS files
image_data_888 = fits.getdata(file_name_888, ext=0)[0][0]
image_data_1367 = fits.getdata(file_name_1367, ext=0)[0][0]
image_data_1419 = fits.getdata(file_name_1419, ext=0)[0][0]

# Plot original images
plt.figure(figsize=(15, 5))
plt.subplot(1, 3, 1)
plt.imshow(image_data_888)
plt.colorbar(cmap='viridis')
plt.title('Original Image (888)')
plt.subplot(1, 3, 2)
plt.imshow(image_data_1367)
plt.colorbar(cmap='viridis')
plt.title('Original Image (1367)')
plt.subplot(1, 3, 3)
plt.imshow(image_data_1419)
plt.colorbar(cmap='viridis')
plt.title('Original Image (1419)')
plt.tight_layout()
plt.show()

#%%
# Plot SNR masked images
snr_cut = 5
# Estimate background noise
bn_888 = estimate_background_noise(image_data_888, 10, 10, width=10, height=10)
bn_1367 = estimate_background_noise(image_data_1367, 10, 10, width=10, height=10)
bn_1419 = estimate_background_noise(image_data_1419, 10, 10, width=10, height=10)

# Calculate SNR images
snr_image_888 = calculate_snr(image_data_888, background_noise=bn_888, neighborhood_radius=2)
snr_image_1367 = calculate_snr(image_data_1367, background_noise=bn_1367, neighborhood_radius=2)
snr_image_1419 = calculate_snr(image_data_1419, background_noise=bn_1419, neighborhood_radius=2)

# Apply mask based on SNR threshold
masked_image_888 = np.where(np.abs(snr_image_888) <= snr_cut, np.nan, image_data_888)
masked_image_1367 = np.where(np.abs(snr_image_1367) <= snr_cut, np.nan, image_data_1367)
masked_image_1419 = np.where(np.abs(snr_image_1419) <= snr_cut, np.nan, image_data_1419)

# Plot masked images
plt.figure(figsize=(15, 5))
plt.subplot(1, 3, 1)
plt.imshow(masked_image_888)
plt.colorbar(cmap='viridis')
plt.title('SNR Maked Image (888)')

plt.subplot(1, 3, 2)
plt.imshow(masked_image_1367)
plt.colorbar(cmap='viridis')
plt.title('SNR Masked Image (1367)')

plt.subplot(1, 3, 3)
plt.imshow(masked_image_1419)
plt.colorbar(cmap='viridis')
plt.title('SNR Masked Image (1419)')

plt.tight_layout()
plt.show()


#%%
# Calculate the spectral index
x = np.log10([888, 1367])
y = np.log10([np.nanmean(masked_image_888), np.nanmean(masked_image_1367)])

slope, intercept, r_value, p_value, std_err = linregress(x, y)
spectral_index = slope

# Calculate the spectral index error
flux_error_888 = np.nanstd(masked_image_888)
flux_error_1367 = np.nanstd(masked_image_1367)
flux_error_1419 = np.nanstd(masked_image_1419)

# Covariance matrix
cov_matrix = np.cov(np.vstack([masked_image_888.flatten(), masked_image_1367.flatten(), masked_image_1419.flatten()]))

# Spectral index uncertainty
spectral_index_error = np.sqrt(cov_matrix[0, 0] * slope ** 2 + cov_matrix[1, 1] + 2 * cov_matrix[0, 1] * slope) / (np.log10(1419) - np.log10(888))

# Plot the masked spectral index map and the best-fit line
plt.imshow(masked_image_888, cmap='jet')
plt.colorbar(label='Intensity')
plt.title('Masked Image (888)')
plt.show()

plt.imshow(masked_image_1367, cmap='jet')
plt.colorbar(label='Intensity')
plt.title('Masked Image (1367)')
plt.show()

plt.imshow(masked_image_1419, cmap='jet')
plt.colorbar(label='Intensity')
plt.title('Masked Image (1419)')
plt.show()

plt.plot(x, y, 'bo', label='Data Points')
plt.plot(x, slope * x + intercept, 'r-', label='Best Fit Line')
plt.xlabel('log(Frequency)')
plt.ylabel('log(Flux)')
plt.title('Spectral Index')
plt.legend()
plt.show()

print('Spectral Index: {:.3f}'.format(spectral_index))
print('Spectral Index Error: {:.3f}'.format(spectral_index_error))


# Print the best-fit line parameters and the spectral index error
# print("Best Fit Line Parameters:")
# print("Slope:", slope)
# print("Intercept:", intercept)
# print("Standard Error:", std_err)
# print("Spectral Index Error:", spectral_index_error)


# %%
header = fits.Header()
header['SIMPLE'] = True
header['BITPIX'] = -64
header['NAXIS'] = 2
header['NAXIS1'] = spectral_index_masked.shape[1]
header['NAXIS2'] = spectral_index_masked.shape[0]
header['BUNIT'] = 'Spectral Index'
header['CTYPE1'] = 'RA---SIN'
header['CTYPE2'] = 'DEC--SIN'
header['CRPIX1'] = 1067
header['CRPIX2'] = 819
header['CRVAL1'] = 82.5
header['CRVAL2'] = -68.6714
header['CDELT1'] = -2.5 / 3600  # Convert from arcsec to degrees
header['CDELT2'] = 2.5 / 3600  # Convert from arcsec to degrees
header['CUNIT1'] = 'deg'
header['CUNIT2'] = 'deg'
header['EQUINOX'] = 2000

hdu = fits.PrimaryHDU(data=spectral_index_masked, header=header)

output_file = 'spectral_index_new_2.fits'
hdu.writeto(output_file, overwrite=True)
# %%
