import numpy as np
from astropy.io import fits
import glob
from astropy.io.fits import getval
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd


list = glob.glob('/Volumes/VINCE/dwarfs/Gonly/*.fits*')
galaxies = pd.read_csv('galaxies.txt', sep=r'\s+')
catalog  = SkyCoord(galaxies.RA * u.degree, galaxies.DEC * u.degree)

df = pd.DataFrame(columns = ['ID', 'RA', 'DEC', 'band', 'file'], index = range(0,len(list)))

for j , file in enumerate(list):
    
    hdu = fits.open(file)
    data = hdu[0].data
    header = hdu[0].header

    # fix filter keyword
    
    filt = header['FILTER']
    if filt == '':
	   filt = file.split('.')[1].lower() # some files are missing the filter keyword
    else:
	   filt = filt[0].lower()

    # get the number of combined images
    history = np.array(header['HISTORY'])
    mask = np.array([history[i].startswith('=') for i in range(len(history))])
    ncomb = len(history[mask])

    w = WCS(file)
    ra, dec = w.wcs_pix2world(648/2., 648/2., 1)

    single = SkyCoord(ra=ra*u.degree, dec=dec*u.degree) 
    idx, d2d, _ = single.match_to_catalog_sky(catalog)

    galaxy_id = galaxies.iloc[idx.item()].ID

    df.loc[j] = pd.Series({'ID':galaxy_id, 'RA':ra, 'DEC':dec, 'band':filt, 'file': file})

    header.set('EXPTIME', 1, 'Exposure time')
    header.set('NCOMBINE', ncomb , 'Number of combined images')
    header.set('GALAXY', galaxy_id , 'ID from Davies+15')
    header.set('GALAXY_RA', ra , 'RA from Davies+15')
    header.set('GALAXY_DEC', dec , 'DEC from Davies+15')

    hdu.writeto('/Volumes/VINCE/dwarfs/final/{}_{}.fits'.format(galaxy_id,j))

    print j, galaxy_id, ra, dec

dup = df[df['ID'].duplicated(keep=False)].sort('ID')

