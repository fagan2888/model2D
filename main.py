from astropy.io import fits
import glob
from astropy.wcs import WCS, utils
import ntpath
from model2D import model2D
import subprocess
import os
from tqdm import *
import pandas as pd

working_dir = os.getcwd() + '/' # Needed to copy psf files

def cleanfiles(listfiles, clean = True):
	ext = ['out.fits','ldac','galfit.txt','mask','.psf', 'psf_', '_sub']
	if clean:
		for f in listfiles:
			for e in ext:
				if e in f:
					os.remove(f)

## Delete files from previous modelling
listfiles = glob.glob('/Volumes/VINCE/dwarfs/combined/*')
cleanfiles(listfiles, clean = True)

# Put actual data files into a list
#listfiles = glob.glob('/Volumes/VINCE/dwarfs/combined_LSBVCC/*.fits*')
#listfiles = pd.read_csv('listfiles.csv').values.flatten()

for f in tqdm(listfiles):
	hdu = fits.open(f)                 # Load the multi-extension fits file

	for i in range(0,len(hdu)):		
		band = hdu[i].name

		if band == 'G':                # Look only for G-band observations

			header = hdu[i].header
			w = WCS(hdu[i].header)

			# Define class parameters
			pscale = str(utils.proj_plane_pixel_scales(w)[0] * 3600)[0:5]
			gain = str(header['GAIN'])
			seeing = str(header['FINALIQ'])
			sizeX = str(header['NAXIS2'])
			sizeY = str(header['NAXIS1'])
			saturate = str(header['SATURATE'])
	
			output = f.split('.fits')[0] + '_{}'.format(band)
			image = f + '[{}]'.format(i)

			# Inizialite the class
			c = model2D(image, output, gain, pscale, seeing, 
				        saturate, sizeX, sizeY, model2D.params, isofactor = 3.)
		
			# Create median subtracted image
			_ = c.get_median()

			# Extract sources from median subtracted image
			astrotable = c.run_sex(output + '_sub.fits')
			
			# Create bad-pixel mask
			_ = c.make_mask(astrotable, flux_radius = 20, fwhm_image = 20)
			
			# Extract sources from original image and create ldac file
			astrotable = c.run_sex(image)
			
			# Create psf image using ldac file created above
			_ = c.run_psfex()

			## Copy psf_XXX.fits files from working directory to data directory
			base, _ = os.path.split(f)
			src = working_dir + 'psf_{}'.format(ntpath.basename(output)) + '.fits'
			dst = base + '/' + 'psf_{}'.format(ntpath.basename(output)) + '.fits'
			os.rename(src,dst)
			##

			# Create galfit input file
			_ = c.make_galfit(f, i)

	hdu.close()

# ------------------------------------------------------


listfiles = glob.glob('/Volumes/VINCE/dwarfs/combined/*galfit*')
for f in l2:
	subprocess.call(['/Volumes/VINCE/dwarfs/analysis/./galfit', f ])







