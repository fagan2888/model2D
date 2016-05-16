import numpy as np
from astropy.io import fits
import glob
from astropy.io.fits import getval, getheader, getdata
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd

datain = '/Volumes/VINCE/dwarfs/data_EVCC/*.fits*'
catin = '/Volumes/VINCE/dwarfs/EVCC/EVCC_total.csv'

files = glob.glob(datain)                 # Directory with the data
galaxies = pd.read_csv(catin, sep=',')    # Object catalogue
catalog  = SkyCoord(galaxies.RA * u.degree, galaxies.DEC * u.degree)  

df = pd.DataFrame(files,columns=['files']) # Copy files to a dataframe


df['suff'] = df['files'].str.split('__').str[1] # Create a column with unique 
												# identification code: the number following 
												# NGVS-X-X.X.0__ 

grouped = df.groupby('suff') # group by unique identification code

for name, group in grouped:
    
	single_files = group['files']
	num = len(single_files)

	hdulist = fits.HDUList()
	for i in range(0,num):
			f = single_files.values[i]

			hdu = fits.open(f)
			header = hdu[0].header

			filt = f.split('.')[1].lower() 

			history = np.array(header['HISTORY'])
			mask = np.array([history[i].startswith('=') for i in range(len(history))])
			ncomb = len(history[mask])

			w = WCS(f)
			size = header['NAXIS1']
			ra, dec = w.wcs_pix2world(size/2., size/2., 1)

			single = SkyCoord(ra = ra*u.degree, dec = dec*u.degree) 
			idx, d2d, _ = single.match_to_catalog_sky(catalog)

			galaxy_id = galaxies.iloc[idx.item()].ID
			GALAXY_RA = galaxies.iloc[idx.item()].RA
			GALAXY_DEC = galaxies.iloc[idx.item()].DEC

			header.set('EXPTIME', 1, 'Exposure time')
			header.set('NCOMBINE', ncomb , 'Number of combined images')
			header.set('GALAXY', galaxy_id , 'ID from Davies+15')
			header.set('GALAXY_RA', GALAXY_RA , 'RA from Davies+15')
			header.set('GALAXY_DEC', GALAXY_DEC , 'DEC from Davies+15')

			insert = fits.ImageHDU(data = hdu[0].data, header=hdu[0].header, name = filt)
			hdu.close()

			hdulist.append(insert)

	hdulist.writeto('/Volumes/VINCE/dwarfs/combined_EVCC/{}'.format(name), clobber = True)




