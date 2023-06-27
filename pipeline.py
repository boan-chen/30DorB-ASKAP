#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import generic_filter
import astropy.io.fits as fits

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
def calculate_snr(image_data, background_noise=0, neighborhood_radius=3):
    signal = image_data - background_noise
    noise = local_std(image_data, neighborhood_radius)
    snr = np.divide(signal, noise, out=np.zeros_like(signal), where=noise != 0)
    return snr

#%%
# Directory and file names
# directory = '/Users/Brian/Documents/2022 Fall Semester/Insterstellar Medium/ProjectData/'
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

#%% Plot SNR masked images
snr_cut = 5
# Estimate background noise
bn_888 = estimate_background_noise(image_data_888, 5, 5, width=10, height=10)
bn_1367 = estimate_background_noise(image_data_1367, 5, 5, width=10, height=10)
bn_1419 = estimate_background_noise(image_data_1419, 5, 5, width=10, height=10)

# Calculate SNR images
snr_image_888 = calculate_snr(image_data_888, background_noise=bn_888)
snr_image_1367 = calculate_snr(image_data_1367, background_noise=bn_1367)
snr_image_1419 = calculate_snr(image_data_1419, background_noise=bn_1419)

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

#%% Plot masked spectral index
spectral_index_snr_cut = 0
# Calculate flux values and errors
flux_888 = masked_image_888
flux_1367 = masked_image_1367
flux_1419 = masked_image_1419

# Calculate flux errors using SNR
flux_error_888 = np.abs(flux_888) / snr_image_888
flux_error_1367 = np.abs(flux_1367) / snr_image_1367
flux_error_1419 = np.abs(flux_1419) / snr_image_1419

# Calculate the spectral index
spectral_index = np.log10(flux_888 / flux_1367) / np.log10(888 / 1367)

# Calculate the spectral index errors (error propagation)
flux_ratio_error =  np.sqrt((flux_error_888 / flux_888) ** 2 + (flux_error_1367 / flux_1367) ** 2)


# spectral_index_error = 0.434 * flux_ratio_error * np.abs(flux_888 / flux_1367) / np.log10(888 / 1367)
# spectral_index_snr = np.abs(spectral_index) / spectral_index_error

# Mask the spectral index values with absolute error <= 3
spectral_index_masked = np.where(np.abs(1/flux_ratio_error) <= spectral_index_snr_cut, np.nan, spectral_index)

plt.imshow(spectral_index_masked, cmap='jet', vmin=-2, vmax=2)
plt.colorbar(label='Spectral Index')
plt.show()
#%% Plot histogram

# Flatten the spectral index masked array and filter out NaN values
spectral_index_masked_flat = spectral_index_masked.flatten()
spectral_index_masked_flat = spectral_index_masked_flat[~np.isnan(spectral_index_masked_flat)]

# Compute the histogram
hist, bins = np.histogram(spectral_index_masked_flat, bins=40, density=True)

# Filter histogram bins with counts > threshold (0.015 in this case)
threshold = 0.05
filtered_bins = bins[:-1][hist > threshold]
filtered_hist = hist[hist > threshold]

# Calculate the bin width
bin_width = np.diff(bins)[0]

# Plot the filtered and normalized histogram
plt.bar(filtered_bins, filtered_hist, width=bin_width, color='blue', alpha=0.7)
plt.xlabel('Spectral Index')
plt.ylabel('Normalized Frequency')
plt.title(f'Normalized Histogram of Spectral Index (Counts > {threshold})')
plt.grid(True)
plt.show()

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

output_file = 'spectral_index_new.fits'
hdu.writeto(output_file, overwrite=True)
# %%
