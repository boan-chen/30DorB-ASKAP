#%%test

import matplotlib.pyplot as plt
from astropy.io import fits

# Open the FITS file
fits_file = 'spectral_index.fits'
hdul = fits.open(fits_file)

# Read the spectral index data and header
spectral_index_data = hdul[0].data
header = hdul[0].header

# Close the FITS file
hdul.close()

# Plot the spectral index data
plt.imshow(spectral_index_data, cmap='jet', vmin=-2, vmax=2)
plt.colorbar(label='Spectral Index')
plt.title('Spectral Index Map')
plt.show()

# %%
