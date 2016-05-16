from astropy.io.fits import getdata
import glob
import shutil

listfiles = glob.glob('/Volumes/VINCE/dwarfs/data_EVCC/*.fits*')

for files in listfiles:

	data = getdata(files)
	zeros = data[data == 0]
	if len(zeros) > 300000:
		print len(zeros)
		shutil.move(files, '/Volumes/VINCE/dwarfs/zeros_EVCC/')
