from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Retrieving cube data from fits file of first exposure
hdul1 = fits.open('Spectra/a01800_/OS.20250720.39412.73.lev1.fits.gz')
header = hdul1[0].header
data_cube1 = hdul1[0].data
data_cube1 = np.transpose(data_cube1, axes=(2, 1, 0))  # Orients for easy viewing

median_cube1 = np.flipud(np.median(data_cube1, axis=0))  # Median intensity 2D image

# Retrieving cube data from fits file of second exposure
hdul2 = fits.open('Spectra/a01800_/OS.20250720.40383.26.lev1.fits.gz')
data_cube2 = hdul2[0].data
data_cube2 = np.transpose(data_cube2, axes=(2, 1, 0))  # Orients for easy viewing

median_cube2 = np.flipud(np.median(data_cube2, axis=0))  # Median intensity 2D image

print(data_cube1.shape)
print(data_cube2.shape)

assert data_cube1.shape[0] == data_cube2.shape[0], \
    'Incompatible wavelength range'

assert data_cube1.shape[2] == data_cube2.shape[2], \
    'Incompatible spaitial width'

# Identifying dimensions needed to accomodate overlapping images
lambda_range = data_cube1.shape[0]
col_max = data_cube1.shape[2]
row_max = data_cube1.shape[1] + data_cube2.shape[1]
# Creating empty array for combining images
empty_array = np.zeros((lambda_range, row_max, col_max))

#TODO: Retrieve offsets from fits headers
offset = 40  # Temporarily hardcoded in pixels
offset_half = int(offset / 2)


empty_array[:, 32+offset_half:96+offset_half, :] = data_cube1
empty_array[:, 32-offset_half:96-offset_half, :] = empty_array[:, 32-offset_half:96-offset_half, :] + -data_cube2
blah = np.flipud(np.median(empty_array, axis=0))
plt.imshow(blah, cmap='gray')
plt.show()
plt.clf()

new_hdu = fits.PrimaryHDU(empty_array, header=header)
new_hdu.writeto('Spectra/a01800_/combined.fits', overwrite=True)

plt.imshow(median_cube1, cmap='gray')
plt.title('Prime Target 1 (a01800_)')
plt.show()
plt.clf()

plt.imshow(median_cube2, cmap='gray')
plt.title('Prime Target 2 (a01800_)')
plt.show()
plt.clf()
