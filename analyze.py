import numpy as np
from astropy.io import fits
import glob
from astropy.io.fits import getval, getdata
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, Normalize
import matplotlib.gridspec as gridspec
from tqdm import *
import re
from scipy import ndimage
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

def sector(header):
    '''
    It extracts the 'image region to fit' in the galfit input file (row H)
    
    INPUT 
    'header' is the second extension of the galfit output file. 

    OUTPUT
    It returns a tuple with (xmin, xmax, ymin, ymax)
    '''
    FITSECT = header['FITSECT']
    sect = re.split('[][:,]',FITSECT)
    x1,x2,y1,y2 = int(sect[1]), int(sect[2]), int(sect[3]), int(sect[4])
    return x1, x2, y1, y2


def plot_model(galaxy, sect = [x1, x2, y1, y2], directory = '/Volumes/VINCE/dwarfs/combined_VCC/figures/'):
    '''
    A wrapper to plot galfit results, bad pixel mask and other info
    
    INPUT
    'galaxy': Single row of the dataframe produced by the procedure the for loop
              below
    'sect'  : Image sector produced by sector(header). It is the same for all 
              galaxies. Do not need to call sector(header) every time.
    'directory': Where to save the results. Directory needs to be created
    '''
    plt.ioff()
    fig, axarr = plt.subplots(2,3)
    #fig.suptitle('{}'.format(f))

    hdu = fits.open(galaxy['MODEL'])
    image = ndimage.gaussian_filter(hdu[1].data, 1)
    axarr[0, 0].imshow(image, cmap='gray', norm=LogNorm(), vmin=1, vmax = 50)
    axarr[0, 0].set_title('Image')

    model = hdu[2].data
    axarr[0, 1].imshow(model, cmap='Blues', norm=LogNorm(), vmin=0.01, vmax = 6)
    axarr[0, 1].set_title('Model')

    residuals = ndimage.gaussian_filter(hdu[3].data, 1)
    axarr[1, 0].imshow(residuals, cmap='gray', norm=LogNorm(), vmin=1, vmax = 50)
    axarr[1, 0].set_title('Residuals')

    x1, x2, y1, y2 = sect[0], sect[1], sect[2], sect[3]
    themask = themask = getdata(galaxy['MASK'])[x1:x2,y1:y2]
    axarr[1, 1].imshow(themask, cmap='gray', vmin=0, vmax = 1)
    axarr[1, 1].set_title('Mask')
    
    axarr[0, 2].text(0.2, 1.0,r'Galaxy: VCC{}'.format(galaxy['ID']), va="center", ha="left")
    axarr[0, 2].text(0.2, 0.9,r'mtot $=$ {}'.format(galaxy['mtot']), va="center", ha="left")
    axarr[0, 2].text(0.2, 0.8,r'Re $=$ {} pc'.format(galaxy['Re']), va="center", ha="left")
    axarr[0, 2].text(0.2, 0.7,r'n $=$ {}'.format(galaxy['n']), va="center", ha="left")
    axarr[0, 2].text(0.2, 0.6,r'PA $=$ {}'.format(galaxy['PA']), va="center", ha="left")
    axarr[0, 2].text(0.2, 0.5,r'chi2nu $=$ {}'.format(galaxy['chi2nu']), va="center", ha="left")
    
    axarr[0, 2].set_aspect('equal')
    axarr[0, 2].axis('off')

    axarr[1, 2].set_aspect('equal')
    axarr[1, 2].axis('off')

    fig.subplots_adjust(hspace=0., wspace = 0.)
    plt.setp([a.get_xticklabels() for a in axarr.flatten()], visible=False);
    plt.setp([a.get_yticklabels() for a in axarr.flatten()], visible=False);

    fig.savefig(directory + '{}.png'.format(f.split('/')[-1]), dpi= 200)
    plt.close(fig)

################
################

dis = 16.5
kpc = (dis * 1000/206265.)
pc = (dis * 1000/206265.) * 1000.

galaxies = pd.read_csv('/Volumes/VINCE/dwarfs/VCC/VCC_clean.csv', sep=',')
l1 = glob.glob('/Volumes/VINCE/dwarfs/combined_VCC/*out.fits*')
l2 = glob.glob('/Volumes/VINCE/dwarfs/combined_LSBVCC/*out.fits*')

listfiles = l1 + l2
listfiles = listfiles[0:10]

columns = ['ID', 'XC', 'YC', 'mtot', 'Re', 
            'n', 'ba','PA', 'chi2nu', 'FINALIQ','MODEL','MASK']

df = pd.DataFrame(columns = columns, index = range(0,len(listfiles)))

for j , f in tqdm(enumerate(listfiles)):
    print f
    hdu = fits.open(f)

    header1 = hdu[1].header
    header2 = hdu[2].header

    #w = WCS(hdu[1].header)
    #RA, DEC = w.all_pix2world(30, 40, 0)

    hdu.close()

    x1, x2, y1, y2 = sector(header2)
    mask_filename = f.split('_out.fits')[0] + '_mask.fits'
    themask = getdata(mask_filename)[x1:x2,y1:y2]
  
    df.loc[j] = pd.Series({
        'ID':header1['GALAXY'], 
        'XC': header2['2_XC'].split(' ')[0], 
        'YC': header2['2_YC'].split(' ')[0], 
        'mtot': header2['2_MAG'].split(' ')[0], 
        'Re': header2['2_RE'].split(' ')[0], 
        'n': header2['2_N'].split(' ')[0], 
        'ba': header2['2_AR'].split(' ')[0], 
        'PA': header2['2_PA'].split(' ')[0], 
        'chi2nu': header2['CHI2NU'],
        'FINALIQ': header1['FINALIQ'],
        'MODEL': f,
        'MASK': mask_filename})

    #######################################

for i in df.index:
    plot_model(df.iloc[i])

    #######################################

df.mtot = df.mtot.astype(float) 
df.re = df.Re.astype(float) 
df.n = df.n.astype(float)
df.ba = df.ba.astype(float)
df.PA = df.PA.astype(float)
df.chi2nu = df.chi2nu.astype(float)

merged = pd.merge(galaxies, df, on = 'ID')

merged.to_csv('result_VCC.csv', index = False)

VCC = pd.read_csv('/Volumes/VINCE/dwarfs/analysis/result_VCC.csv')
LSBVCC = pd.read_csv('/Volumes/VINCE/dwarfs/analysis/result_LSBVCC.csv')


cm = plt.cm.get_cmap('RdYlBu')
fig = plt.figure()
plt.gca().invert_xaxis()
plt.yscale('log')
plt.scatter(df.mtot, df.n, c=df.chi2nu, s=35, cmap=cm)

plt.scatter(VCC.mtot, VCC.Re, s=35, c = 'red')
