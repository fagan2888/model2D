from astropy.io import fits
import glob
from pyraf import iraf

listfiles = glob.glob('/Volumes/VINCE/dwarfs/zeros_EVCC/*.fits*')
iraf.images()

for files in listfiles:

	iraf.imcopy(input = files + '[2:968,2:968]', output = files )
	
	hdu = fits.open(files)

	band = hdu[0].name
	header = hdu[0].header

	sizeX = str(header['NAXIS2'])
	sizeY = str(header['NAXIS1'])

	hdu.close()