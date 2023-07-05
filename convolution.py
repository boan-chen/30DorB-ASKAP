#%%
import numpy as np
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian2DKernel
from astropy.io import fits

# Load the original image
file_dictionary = {'888': 'data/signal/concolved_887.491MHz_18arcsec_30dorB_signal.fits', \
                   '1367': 'data/signal/concolved_1367.49MHz_18arcsec_30dorB_signal.fits', \
                   '4800': 'data/signal/LMC4.8_ATCA_signal.fits', \
                   '8600': 'data/signal/LMC8.6_ATCA_signal.fits', \
                   '888_noise': 'data/background/concolved_887.491MHz_18arcsec_30dorB_background.fits', \
                   '1367_noise': 'data/background/concolved_1367.49MHz_18arcsec_30dorB_background.fits', \
                   '4800_noise': 'data/background/LMC4.8_ATCA_background.fits', \
                   '8600_noise': 'data/background/LMC8.6_ATCA_background.fits'}

original_image_file = file_dictionary['888']
original_image_data, original_image_header = fits.getdata(original_image_file, header=True)

if sum(np.shape(np.shape(original_image_data))) == 2:
    cube = False
else:
    cube = True
    
# Get the pixel scale in arcseconds
pixel_scale = np.abs(original_image_header['CDELT2']) * 3600  # Assuming the pixel increment is the same in both axes

# Convert the restoring beam size from arcseconds to pixels
restoring_beam_arcsec = 35  # Desired restoring beam size in arcseconds
restoring_beam_pixels = restoring_beam_arcsec / pixel_scale

# Create the Gaussian kernel
gaussian_kernel = Gaussian2DKernel(restoring_beam_pixels)

# Convolve the image with the Gaussian kernel
if cube == True:
    convolved_image_data_raw = convolve_fft(original_image_data[0][0], gaussian_kernel._array)
    convolved_image_data = original_image_data.copy()
    convolved_image_data[0][0] = convolved_image_data_raw
else:
    convolved_image_data = convolve_fft(original_image_data, gaussian_kernel._array)
    
print(np.shape(gaussian_kernel._array))
print(np.shape(convolved_image_data))
# Save the convolved image to a FITS file
# convolved_image_file = 'data/convolved_signal/convolved_888_signal.fits'
# fits.writeto(convolved_image_file, convolved_image_data, header=original_image_header, overwrite=True)

# %%
