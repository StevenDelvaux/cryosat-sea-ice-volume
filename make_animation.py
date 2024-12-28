import contextlib
from datetime import datetime, timedelta
from PIL import Image

def makeAnimation(enddate, frames, animationFileName, getFileNameFromDate, missingDates = [], endpause = 5):
	date = datetime(enddate.year, enddate.month, enddate.day)
	
	filenames = []
	counter = 0
	while counter < frames:
		print('plotting date: ',date)
		if date in missingDates:
			date = date - timedelta(days = 1)
			print('missing date: ', date)
			continue
		localfilename =  getFileNameFromDate(date)
		filenames.append(localfilename)
		counter += 1
		date = date - timedelta(days = 1)
	filenames.reverse()
	
	lastfilename = filenames[-1]
	for k in range(endpause):
		filenames.append(lastfilename)
	with contextlib.ExitStack() as stack:
		imgs = (stack.enter_context(Image.open(f)) for f in filenames)
		img = next(imgs)
		# https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
		startdate = enddate - timedelta(days = frames-1)
		img.save(fp=animationFileName, format='GIF', append_images=imgs, save_all=True, duration=500, loop=0) #, quality=25, optimize=True, compress_level=9)
		compress_string = "magick mogrify -layers Optimize -fuzz 7% " + animationFileName
