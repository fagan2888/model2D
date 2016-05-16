__version__ = '0.1'

import numpy as np 
import pandas as pd
from astropy.io import fits
import astromatic_wrapper as aw

from astropy.wcs import WCS, utils
from astropy.io.fits import getheader
import ntpath
from pyraf import iraf
import os
iraf.imfilter()

class model2D:

	params = ['X_IMAGE', 'Y_IMAGE', 'ALPHA_J2000','DELTA_J2000','CXX_IMAGE','CYY_IMAGE',
	          'CXY_IMAGE','ISOAREA_IMAGE','ISOAREAF_IMAGE','FWHM_IMAGE','CLASS_STAR',
	          'VIGNET(35,35)','FLUX_RADIUS','FLUX_APER','FLUXERR_APER','ELONGATION',
	          'FLAGS','SNR_WIN', 'MAG_AUTO']
 
	def __init__(self, image, output, gain, pscale, seeing, saturate, sizeX, sizeY, params,
				 isofactor = 2.0, fluxT = 'median', radiusT = 20, verbose = True):

		self.image = image
		self.output = output
		self.gain = gain
		self.pscale = pscale
		self.seeing = seeing
		self.saturate = saturate

		self.sizeX = sizeX
		self.sizeY = sizeY
		self.params = params

		self.isofactor = isofactor
		self.fluxT = fluxT
		self.radiusT = radiusT

		self.verbose = verbose

	def header(self):
		return getheader(self.image)

	def run_sex(self, input_image, verbose = None):
		'''
		Run SExtractor on an astronomical image. 
		The output is the filename with extension .ldac

		params: SExtractor parameters to be computed.
		The full set of parameters available can be printed on screen by typing 'sex -dp' from shell
		'''

		files = {'image' : input_image}
		kwargs = {'code': 'SExtractor'}

		kwargs['config_file'] = 'model2D.sex'

		kwargs['config'] = {'CATALOG_NAME':'ldac.fits'}

		ldacfile = '{}.ldac'.format(self.output)

		kwargs['config']['CATALOG_NAME'] = ldacfile
		kwargs['config']['CATALOG_TYPE'] = 'FITS_LDAC'
		kwargs['config']['FILTER'] = 'Y'
		kwargs['config']['FILTER_NAME'] = 'model2D.conv'
		kwargs['config']['WEIGHT_TYPE'] = 'NONE'
		kwargs['config']['PHOT_APERTURES'] = '3,5,7,10,13,15,18,20'
		kwargs['config']['PHOT_AUTOPARAMS'] = '2.5,3.5'
		kwargs['config']['SATUR_LEVEL'] = self.saturate
		kwargs['config']['MAG_ZEROPOINT'] = '30'
		kwargs['config']['PIXEL_SCALE'] = self.pscale
		kwargs['config']['SEEING_FWHM'] = self.seeing
		kwargs['config']['STARNNW_NAME'] = 'model2D.nnw'
		kwargs['config']['CHECKIMAGE_TYPE'] = 'NONE'

		kwargs['temp_path'] = '.'
		kwargs['params'] = self.params
		
		if verbose == None:
			verbose = self.verbose
			kwargs['store_output'] = verbose

			sextractor = aw.api.Astromatic(**kwargs)
			sextractor.run(files['image'], store_output = verbose)

		astrotable = aw.utils.ldac.get_table_from_ldac(ldacfile)
	
		#if dataframe:
		#	astrotable = pd.DataFrame.from_records(astrotable, columns = astrotable.columns)

		return astrotable

	def get_median(self):
		'''
		Calculates the median subtracted image used to compute the pixel mask
		'''
		iraf.median(input = self.image, output = 'tmp_med.fits', 
					xwindow = 40, ywindow = 40, verbose ='No')

		iraf.imarith(operand1 = self.image, operand2 = 'tmp_med.fits[0]', 
					 op = '-', result = self.output + '_sub.fits')
		
		os.remove('tmp_med.fits')


	def run_psfex(self, verbose = None, copy = True):
		'''
		Run PSFex on an astronomical image. 
		The output is the filename with extension .fits.ldac

		ldacfile: The .ldac filename produced by sextractor.
		'''

		ldacfile = '{}.ldac'.format(self.output)

		kwargs = {'code': 'PSFEx'}
		kwargs['config_file'] = 'model2D.psfex'

		kwargs['config'] = {'BASIS_TYPE':'PIXEL_AUTO'}
		kwargs['config']['PSFVAR_GROUPS'] = '1,1'
		kwargs['config']['PSFVAR_DEGREES'] = '0'
		kwargs['config']['PSFVAR_NSNAP'] = '1'
		kwargs['config']['SAMPLE_FWHMRANGE'] = '2.0,8.0'
		kwargs['config']['CHECKPLOT_DEV'] = 'NULL'

		kwargs['config']['CHECKIMAGE_TYPE'] = 'SNAPSHOTS_IMRES'
		kwargs['config']['CHECKIMAGE_NAME'] = 'psf'

		if verbose == None:
			verbose = self.verbose
			kwargs['store_output'] = verbose
		
			psfex = aw.api.Astromatic(**kwargs)
			psfex.run(ldacfile, store_output = verbose)

			header = getheader('{}.psf'.format(self.output), 'PSF_DATA')

		if copy == True:
			pass


	def make_mask(self, astrotable, inradius = 30, flux_radius = 20, fwhm_image = 15):
		'''
		Make a mask file for galfit. 
		'''
		X = np.linspace(0, int(self.sizeX), int(self.sizeX))
		Y = np.linspace(0, int(self.sizeY), int(self.sizeY))[:,None]
		
		astrotable['R'] = (np.sqrt((astrotable['X_IMAGE'] - int(self.sizeX)/2.)**2 + 
							(astrotable['Y_IMAGE'] - int(self.sizeY)/2.)**2))

		bolmask = ((astrotable['R'] < inradius) & (astrotable['FWHM_IMAGE'] > fwhm_image))
		astrotable = astrotable[~bolmask]

		if self.fluxT == 'median':
			cutoff = np.median(astrotable['FLUX_APER'])

		clean = astrotable[astrotable['FLUX_APER'] > cutoff]
		
		X_IMAGE = clean['X_IMAGE'] 
		Y_IMAGE = clean['Y_IMAGE']
		CXX_IMAGE = clean['CXX_IMAGE']
		CYY_IMAGE = clean['CYY_IMAGE']
		CXY_IMAGE = clean['CXY_IMAGE']
		ISOAREA_IMAGE = clean['ISOAREA_IMAGE']

		ell = np.zeros((len(X),len(Y)))

		for i in range(0,len(X_IMAGE)):

			single = ((CXX_IMAGE[i] * (X - X_IMAGE[i])**2) + (CYY_IMAGE[i] * (Y - Y_IMAGE[i])**2) + 
			    		(CXY_IMAGE[i] * (X - X_IMAGE[i])*(Y - Y_IMAGE[i]))) <= ISOAREA_IMAGE[i]**(1/self.isofactor)

			single = single.astype(int)

			ell += single

		hdumask = fits.PrimaryHDU(ell)
		hdulist = fits.HDUList([hdumask])

		hdulist.writeto('{}_mask.fits'.format(self.output), clobber = True)
		

	def make_galfit(self, input_file, hdu):

		fout = open('{}_galfit.txt'.format(self.output),'w')

		x0 = str(int(self.sizeX)/2.)
		y0 = str(int(self.sizeY)/2.)

		base, filename = os.path.split(input_file)
		filename = filename.split('.')[0]
        
		fout.write("""
================================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) """+base+"""/"""+filename+""".fits["""+str(hdu)+"""]    # Input data image (FITS file)
B) """+base+"""/"""+filename+"""_G_out.fits      # Output data image block
C) none                  	     # Sigma image name (made from data if blank or "none") 
D) """+base+"""/psf_"""+filename+"""_G.fits      # Input PSF image and (optional) diffusion kernel
E) 1                          # PSF fine sampling factor relative to data 
F) """+base+"""/"""+filename+"""_G_mask.fits     # Bad pixel mask (FITS image or ASCII coord list)
G) none                       # File with parameter constraints (ASCII file) 
H) 250 717 250 717      # Image region to fit (xmin xmax ymin ymax)
I) 100    100                 # Size of the convolution box (x y)
J) 30.000                     # Magnitude photometric zeropoint 
K) """+self.pscale+""" """+self.pscale+"""      # Plate scale (dx dy)   [arcsec per pixel]
O) regular                    # Display type (regular, curses, both)
P) 0                          # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps


# Component number: 1
 0) sky                    #  Component type
 1) 0.0            1       #  Sky background at center of fitting region [ADUs]
 2) 0.000e+00      0       #  dsky/dx (sky gradient in x)     [ADUs/pix]
 3) 0.000e+00      0       #  dsky/dy (sky gradient in y)     [ADUs/pix]
 Z) 1                      #  Skip this model in output image?  (yes=1, no=0)

# Component number: 2
 0) sersic        #  Component type
 1) """+x0+""" """+y0+"""  1 1  #  Position x, y
 3) 19       1    #  Integrated magnitude 
 4) 50       1    #  R_e (effective radius)   [pix]
 5) 0.8      1    #  Sersic index n (de Vaucouleurs n=4) 
 9) 0.7      1    #  Axis ratio (b/a)  
10) 0        1    #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0             #  Skip this model in output image?  (yes=1, no=0)

================================================================================
	""")
		fout.close()
