# coding: latin-1
# Version 2022-12-28
# Gridded sea ice thickness data from  ftp://ftp.awi.de/sea_ice/product/cryosat2_smos/v206/
# Regional mask from  ftp://ftp.awi.de/sea_ice/product/cryosat2/v2p5/nh/l3c_grid/isoweek/

import numpy as np
from netCDF4 import Dataset
from datetime import date, datetime, timedelta
import glob
import csv
import sys
import os
import shutil
import urllib.request
from contextlib import closing
from math import sqrt, sin, cos, pi, floor, isnan
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from decouple import config

import time
import dropbox_client
import upload_to_google_drive
import regional_python_graphs
import make_animation
import get_last_saved_day

thresh = 15.            # Concentration threshold for area/extent (%)
sic_unc = 0.05          # Default concentration uncertainty
grid_spacing_km = 25.   # Default EASE grid spacing
monthNames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
monthLengths = [31,28,31,30,31,30,31,31,30,31,30,31]
ftpFolder = 'ftp://ftp.awi.de/sea_ice/product/cryosat2_smos/v206/nh/'

putOnDropbox = True

class RegionCode: 		# Region codes used in CryoSat auxiliary data
	cab = 1
	beaufort = 2
	chukchi = 3
	ess = 4
	laptev = 5
	kara = 6
	barents = 7
	greenland = 8
	baffin = 9
	stlawrence = 10
	hudson = 11
	caa = 12
	bering = 13
	okhotsk = 14
	
def getClosestRegionCode(row,col):
	"""
	Get the region code for the row and column coordinate from the CryoSat region mask.
    Since the region mask has some missing values (with value -1), we look at neighboring entries if necessary 
	until we find one with a valid region code. 
    """
	
	regionCode = mask[row,col]
	if regionCode != -1:
		return regionCode
	
	radius = 1
	while radius < 10: # look at neighboring entries until we find one with a valid region code
		neighbors = []
		for k in range(radius):
			neighbors.append(mask[row+radius-k,col+k])
			neighbors.append(mask[row-k,col+radius-k])
			neighbors.append(mask[row-radius+k,col-k])
			neighbors.append(mask[row+k,col-radius+k])
			
		neighbors = list(filter(lambda x: x != -1, neighbors)) # remove empty region codes
		if neighbors:
			return max(set(neighbors), key=neighbors.count) # return the region code that appears most often in the list of neighbors
		
		radius += 1 # if no valid region code was found among the closest neighbors, look at neighbors that are one step further away

def rounded(n):
	"""
	Transform a number into a string with 2 decimal digits. 
    """
	return ("{:.2f}".format(n))

def justify(n):
	"""
	Transform a number into a string with 2 decimal digits, right justified. 
    """
	return rounded(n).rjust(8) + ' km³'

def padzeros(n):
	"""
	Left pad a number with zeros. 
    """
	return str(n) if n >= 10 else '0'+str(n)

def getFileName(date):
	startDate = date - timedelta(days = 3)
	endDate = date + timedelta(days = 3)
	return 'W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_' + str(startDate.year) + padzeros(startDate.month) + padzeros(startDate.day)+ '_' + str(endDate.year) + padzeros(endDate.month) + padzeros(endDate.day) + '_' + ('r' if date.year < 2024 or date.year == 2024 and date.month < 6 else 'o') + '_v206_01_l4sit.nc'
	
def getGriddedThickness(date):
	filename = 'data/LATEST/' + getFileName(date)
	if not os.path.isfile(filename):
		filename = download(date)
	f = Dataset(filename, 'r', format="NETCDF4")
	thicknessData = f.variables['analysis_sea_ice_thickness'][:]
	f.close()
	return thicknessData

def download(date):
	"""
	Download Cryosat-SMOS ftp file. 
    """
	
	filename = getFileName(date)
	downloadFilename = filename
	if date.year == 2025 and date.month == 3 and date.day == 25:
		downloadFilename = getFileName(datetime(date.year, date.month, date.day-1))
	downloadFilename = downloadFilename.replace(',','%2C')
	ftpSubfolder = (str(date.year) + "/" + padzeros(date.month)) if (date.year < 2024 or date.year == 2024 and date.month < 6) else 'LATEST'
	fullFtpPath = ftpFolder + ftpSubfolder + "/" + downloadFilename
	localpath = 'data/LATEST/' + filename
	print('downloading file ', fullFtpPath, localpath)
	with closing(urllib.request.urlopen(fullFtpPath)) as r:
		with open(localpath, 'wb') as f:
			shutil.copyfileobj(r, f)
	return localpath	
		
def getLatestDate(csvFileName):
	lastSavedStartDay,lastSavedEndDay = get_last_saved_day.getLastSavedDay(csvFileName)
	lastSavedEndDayString = str(lastSavedEndDay)
	print('inside last saved day', lastSavedStartDay, lastSavedEndDay)
	latestDate = datetime(int(lastSavedEndDayString[0:4]), int(lastSavedEndDayString[4:6]), int(lastSavedEndDayString[6:8]))
	return latestDate
		
def downloadNewFiles():
	dayBeforeYesterday = datetime.today() - timedelta(days = 2)
	csvFileName = 'cryosat-smos-regional-volume.csv'	
	dropbox_client.downloadFromDropbox([csvFileName])
	
	latestDate = getLatestDate(csvFileName)
	
	date = latestDate + timedelta(days = 1)
	outFile = open(csvFileName, 'a', newline='')
	csvFile = csv.writer(outFile)
	found = False
	while date < dayBeforeYesterday:
		print('downloading', date, dayBeforeYesterday)
		filename = ''
		try:
			filename = download(date - timedelta(days = 3))
			found = True
		except:
			print('File not found: ', date)
			break
		date = date + timedelta(days = 1)
		csvFile.writerow(dayvol(filename))
	outFile.close()
	date = date - timedelta(days = 1)
	return date	

def getNsidcLandMask():
	landmask = np.genfromtxt('landmask_nsidc.csv', delimiter=',')	
	n2 = 1
	landmask = landmask[n2:-n2,n2:-n2]
	n = landmask.shape[0]  # n is assumed to be an odd number
	global landmaskcenter
	landmaskcenter = (n-1)/2
	return landmask
	
def insertCryosatDataInNsidcMask(cryosatData, day, year, dummyvalue):
	
	landmask = getNsidcLandMask()
	landmask = landmask*dummyvalue
	landmaskSize = landmask.shape[0]
	counter=1000.0
	rand = (year+day-2000)/200000.0 #random.randrange(1,100000)/10000000.0
	
	for i in range(0,432):
		for j in range(0,432):
			v = cryosatData[0,i,j]
			latitude = lat[i,j]
			longitude = lon[i,j]
			rad = 360*sqrt(2)*sin(pi*(90-latitude)/360)
			y = int(round(landmaskcenter+rad*sin(pi*longitude/180.0)))
			x = int(round(landmaskcenter+rad*cos(pi*longitude/180.0)))
			
			if(x < 0 or y < 0 or x >= landmaskSize or y >= landmaskSize or landmask[x,y] == 0): #or mask[i,j] == -1
				continue				
			if not type(v) is np.float64:
				v = 0
			if(landmask[x,y] == dummyvalue):
				landmask[x,y] = 0
			landmask[x,y] += counter + max(v, rand)

	for x in range(0,landmask.shape[0]):
		for y in range(0,landmask.shape[1]):
			if(landmask[x,y] == 0 or landmask[x,y] == dummyvalue):
				continue
			n = floor(landmask[x,y]/counter)
			landmask[x,y] = (landmask[x,y] - counter*n)/n
			
	return landmask

def getInterpolatedValue(x, y, landmask, dummyvalue):
	radius = 1
	while(radius < 10):
		for k in range(radius):
			other = landmask[x+k,y+radius-k] if x+k < 359 and y+radius-k < 359 else 0
			if(other != dummyvalue and other != 0):
				return other
			other = landmask[x+radius-k,y-k] if x+radius-k < 359 and y-k >= 0 else 0
			if(other != dummyvalue and other != 0):
				return other
			other = landmask[x-k,y-radius+k] if x-k >= 0 and y-radius+k >= 0 else 0
			if(other != dummyvalue and other != 0):
				return other
			other = landmask[x-radius+k,y+k] if x-radius+k >= 0 and y+k < 359 else 0
			if(other != dummyvalue and other != 0):
				return other
		radius += 1
	return dummyvalue

def interpolate(landmask, dummyvalue, anomalyplot):
	mask = landmask.copy()
	
	for x in range(0,landmask.shape[0]):
		for y in range(0,landmask.shape[1]):
			if(mask[x,y] == 0): # land
				if(anomalyplot):
					mask[x,y] = -anomalymax
				continue
			if(mask[x,y] == dummyvalue): # value to be interpolated
				mask[x,y] = getInterpolatedValue(x, y, landmask, dummyvalue)
			if((anomalyplot and abs(mask[x,y]) < 0.001) or ((not anomalyplot) and mask[x,y] < 0.05)): # hide in maps
				mask[x,y] = dummyvalue
			if(anomalyplot and mask[x,y] != dummyvalue and mask[x,y] > anomalymax * 0.99):
				mask[x,y] = anomalymax * 0.99
			if(anomalyplot and mask[x,y] < -anomalymax * 0.99):
				mask[x,y] = -anomalymax * 0.99
			if((not anomalyplot) and mask[x,y] != dummyvalue and mask[x,y] > thicknessmax - 0.06):
				mask[x,y] = thicknessmax - 0.06				
					
	return mask

def plotThickness(landmask,plotTitle,filename,dropboxFilename):
	cdict = {'red':   ((0.0,  0.5, 0.5),
					   (0.001, 0.5, 0.0),
		           	   (0.05, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.15, 0.0, 0.0),					   
					   (0.2, 0.0, 0.2),
				   	   (0.25, 0.2, 0.4),					  
					   (0.3, 0.4, 0.6),
				   	   (0.35, 0.6, 0.8),
					   (0.4, 0.8, 1.0),
					   (0.45, 1.0, 1.0),
				   	   (0.5, 1.0, 1.0),
					   (0.55, 1.0, 1.0),
					   (0.6, 1.0, 0.95),
					   (0.65, 0.95, 0.9),
					   (0.7, 0.9, 0.85),
					   (0.75, 0.85, 0.8),
					   (0.8, 0.8, 0.75),
					   (0.85, 0.75, 0.7),
					   (0.9, 0.7, 0.65),
					   (0.95, 0.65, 0.6),
				       (0.999,  0.6, 1),
                       (1.0,  1, 1)),

         'green':      ((0.0,  0.5, 0.5),
         	           (0.001, 0.5, 0.0),
         	           (0.05, 0.0, 0.1),
					   (0.1, 0.1, 0.25),
					   (0.15, 0.25, 0.4),
					   (0.2, 0.4, 0.55),
					   (0.25, 0.55, 0.7),
					   (0.3, 0.7, 0.85),
					   (0.35, 0.85, 1.0),
					   (0.4, 1.0, 1.0),		
					   (0.45, 1.0, 0.9),				   					   
				   	   (0.5, 0.9, 0.8),
					   (0.55, 0.8, 0.75),
					   (0.6, 0.75, 0.7),
					   (0.65, 0.7, 0.6),
					   (0.7, 0.6, 0.5),
					   (0.75, 0.5, 0.4),
					   (0.8, 0.4, 0.3),
					   (0.85, 0.3, 0.2),
					   (0.9, 0.2, 0.1),
					   (0.95, 0.1, 0.0),
         	           (0.999,  0.0, 1),
                       (1.0,  1, 1)),

         'blue':       ((0.0,  0.5, 0.5),
         	           (0.001, 0.4, 0.4),
         	           (0.05, 0.4, 0.55),
					   (0.1, 0.55, 0.7),
					   (0.15, 0.7, 0.85),
				       (0.2, 0.85, 1.0),	
					   (0.25, 1.0, 0.8),						     
				   	   (0.3, 0.8, 0.6),		
					   (0.35, 0.6, 0.4),		   	   
					   (0.4, 0.4, 0.2),		
					   (0.45, 0.2, 0.0),		   				   
					   (0.5, 0.0, 0.0),
					   (0.55, 0.0, 0.0),
					   (0.6, 0.0, 0.0),
					   (0.7, 0.0, 0.0),
					   (0.8, 0.0, 0.0),
					   (0.9, 0.0, 0.0),
         	           (0.999,  0.0, 0.0),
                       (1.0,  1, 1))}
	kleur = LinearSegmentedColormap('BlueRed1', cdict)
	# next part only serves to get a nicer colormap in the plot
	cbrol = {'red':   ((0.0,  0.0, 0.0),
					   (0.05, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.15, 0.0, 0.0),					   
					   (0.2, 0.0, 0.2),
				   	   (0.25, 0.2, 0.4),					  
					   (0.3, 0.4, 0.6),
				   	   (0.35, 0.6, 0.8),
					   (0.4, 0.8, 1.0),
					   (0.45, 1.0, 1.0),
				   	   (0.5, 1.0, 1.0),
					   (0.55, 1.0, 1.0),
					   (0.6, 1.0, 0.95),
					   (0.65, 0.95, 0.9),
					   (0.7, 0.9, 0.85),
					   (0.75, 0.85, 0.8),
					   (0.8, 0.8, 0.75),
					   (0.85, 0.75, 0.7),
					   (0.9, 0.7, 0.65),
					   (0.95, 0.65, 0.6),		
                       (1.0,  0.6, 0.6)),


			 'green': ((0.0,  0.0, 0.0),
					   (0.05, 0.0, 0.1),
					   (0.1, 0.1, 0.25),
					   (0.15, 0.25, 0.4),
					   (0.2, 0.4, 0.55),
					   (0.25, 0.55, 0.7),
					   (0.3, 0.7, 0.85),
					   (0.35, 0.85, 1.0),
					   (0.4, 1.0, 1.0),		
					   (0.45, 1.0, 0.9),				   					   
				   	   (0.5, 0.9, 0.8),
					   (0.55, 0.8, 0.75),
					   (0.6, 0.75, 0.7),
					   (0.65, 0.7, 0.6),
					   (0.7, 0.6, 0.5),
					   (0.75, 0.5, 0.4),
					   (0.8, 0.4, 0.3),
					   (0.85, 0.3, 0.2),
					   (0.9, 0.2, 0.1),
					   (0.95, 0.1, 0.0),
					   (1.0,  0.0, 0.0)),

			 'blue':  ((0.0, 0.4, 0.4),
					   (0.05, 0.4, 0.55),
					   (0.1, 0.55, 0.7),
					   (0.15, 0.7, 0.85),
				       (0.2, 0.85, 1.0),	
					   (0.25, 1.0, 0.8),						     
				   	   (0.3, 0.8, 0.6),		
					   (0.35, 0.6, 0.4),		   	   
					   (0.4, 0.4, 0.2),		
					   (0.45, 0.2, 0.0),		   				   
					   (0.5, 0.0, 0.0),
					   (0.55, 0.0, 0.0),
					   (0.6, 0.0, 0.0),
					   (0.7, 0.0, 0.0),
					   (0.8, 0.0, 0.0),
					   (0.9, 0.0, 0.0),
					   (1.0,  0.0, 0.0))}
	kleurbrol = LinearSegmentedColormap('BlueRed2', cbrol)
	mask = landmask[30:-70,10:-70]#[50:-90,80:-90]#landmask[85:-100,95:-100]#landmask[30:-70,10:-70]
	n = landmask.shape[0]
	try:
		plt.colorbar().remove()
	except:
		print('error remove color bar thickness')
	plt.clf()
	plt.cla()
	figbrol = plt.imshow(mask, extent=(0,n,0,n), vmin= 0, vmax=thicknessmax,
			   interpolation='nearest', cmap=kleurbrol)

	#plot the relevant map:
	fig2 = plt.imshow(mask, extent=(0,n,0,n), vmin= 0, vmax=thicknessmax,
           interpolation='nearest', cmap=kleur)
	
	plt.title(plotTitle)
	cb = plt.colorbar(figbrol)
	plt.xticks([])
	plt.yticks([])
	cb.set_label("meters")
	plt.savefig(filename)
	if dropboxFilename != '':
		plt.savefig(dropboxFilename + '.png')
	
def plotAnomaly(landmask, plotTitle, filename, dropboxFilename):
	cdict = {'red': ((0.0,  0.4, 0.4),
         	       (0.001, 0.0, 0.0),
         	       #(0.4, 0.8, 0.8),
         	       (0.5, 1.0, 1.0),
         	       (0.999,  0.0, 0.0),
                   (1.0,  1, 1)),

         'green':   ((0.0,  0.4, 0.4),
		           (0.001, 0.0, 0.0),
		           (0.5, 1.0, 1.0),
				   (0.999,  1.0, 1.0),
                   (1.0,  1, 1)),

         'blue':  ((0.0,  0.4, 0.4),
         	       (0.001, 0.4, 0.4),
         	       #(0.4, 1, 0),
         	       (0.5, 1, 0.5),
         	       (0.999,  0.0, 0.0),
                   (1.0,  1, 1))}
	kleur = LinearSegmentedColormap('BlueRed3', cdict)
	# next part only serves to get a nicer colormap in the plot
	cbrol = {'red':   ((0.0,  0.0, 0.0),
					   (0.5, 1, 1),
					   (1.0,  0.0, 0.0)),
			 'green': ((0.0,  0, 0),
					   (0.5, 1, 1),
					   (1.0,  1, 1)),

			 'blue':  ((0.0, 0.4, 0.4),
					   #(0.4, 1, 0),
					   (0.5, 1, 0.5),
					   (1.0,  0.0, 0.0))}
	kleurbrol = LinearSegmentedColormap('BlueRed4', cbrol)
	mask = landmask[30:-70,10:-70]#[50:-90,80:-90]#landmask[85:-100,95:-100]#landmask[30:-70,10:-70]
	n = mask.shape[0]
	m = mask.shape[1]
	try:
		plt.colorbar().remove()
	except:
		print('error remove color bar anomaly')
	plt.clf()
	plt.cla()
	figbrol = plt.imshow(mask, extent=(0,n,0,n), vmin= -anomalymax, vmax=anomalymax,
			   interpolation='nearest', cmap=kleurbrol)

	#plot the relevant map:
	fig2 = plt.imshow(mask, extent=(0,n,0,n), vmin= -anomalymax, vmax=anomalymax,
           interpolation='nearest', cmap=kleur)
	
	plt.title(plotTitle)
	cb = plt.colorbar(figbrol)
	plt.xticks([])
	plt.yticks([])
	cb.set_label("meters")
	#plt.show()
	plt.savefig(filename)
	if dropboxFilename != '':
		plt.savefig(dropboxFilename + '.png')
	
def addMasks(landmask, mask, multiplier, dummyvalue):
	for x in range(0,landmask.shape[0]):
		for y in range(0,landmask.shape[1]):
			if(landmask[x,y] == 0 or landmask[x,y] == dummyvalue):
				continue
			landmask[x,y] = landmask[x,y] + multiplier*mask[x,y]

	return landmask
	
def dayvol(filename) :
	"""
	Calculate regional volume for a daily gridded thickness file. 
    """	
	dates = filename.split('_')
	startstr = dates[5]
	endstr = dates[6]
    	
	f = Dataset(filename, 'r', format="NETCDF4")
	
	# read sea ice concentration, thickness and thickness uncertainty
	
	sic = f.variables['sea_ice_concentration'][:].squeeze()
	sit = f.variables['analysis_sea_ice_thickness'][:]
	sit_unc = f.variables['analysis_sea_ice_thickness_unc'][:].squeeze()	
	
	f.close()

	gg = grid_spacing_km**2   # Area of 25km EASE grid
	gridarea = np.full((432,432), gg)

	# Volume
	per_grid_cell_volume = (sic/100.) * sit * gg
	volume = np.nansum(per_grid_cell_volume) / 1000.0
	per_grid_cell_uncertainty = per_grid_cell_volume * np.sqrt((sic_unc/sic)**2. + (sit_unc/sit)**2.)
	volume_uncertainty = np.nansum(per_grid_cell_uncertainty) / 1000.0

	# Extent and Area
	extent = gridarea[(sic>=thresh)&(sic<=100.)].sum()
	area = gg * 0.01 * sic[(sic>=thresh)&(sic<=100.)].sum()

	# Regional volume
	vtotal, vhudson, vlawrence, vbaffin, vcaa, vbeaufort, vchukchi, vess, vlaptev, vkara, vbarents, vgreenland, vcab, vbering, vokhotsk, vother = (0,)*16
	_,numberOfRows,numberOfColumns = per_grid_cell_volume.shape
	
	for row in range(numberOfRows): 		# numberOfRows = 432
		for col in range(numberOfColumns): 	# numberOfColumns = 432
		
			entry = per_grid_cell_volume[0,row,col]
			if type(entry) is np.float64:
				entry /= 1000.0
				entry = round(entry,3)
				vtotal += entry

				regionCode = getClosestRegionCode(row,col)

				if regionCode == RegionCode.cab:
					vcab += entry
				elif regionCode == RegionCode.beaufort:
					vbeaufort += entry
				elif regionCode == RegionCode.chukchi:
					vchukchi += entry
				elif regionCode == RegionCode.ess:
					vess += entry
				elif regionCode == RegionCode.laptev:
					vlaptev += entry
				elif regionCode == RegionCode.kara:
					vkara += entry
				elif regionCode == RegionCode.barents:
					vbarents += entry
				elif regionCode == RegionCode.greenland:
					vgreenland += entry
				elif regionCode == RegionCode.baffin:
					vbaffin += entry
				elif regionCode == RegionCode.stlawrence:
					vlawrence += entry
				elif regionCode == RegionCode.hudson:
					vhudson += entry
				elif regionCode == RegionCode.caa:
					vcaa += entry
				elif regionCode == RegionCode.bering:
					vbering += entry
				elif regionCode == RegionCode.okhotsk:
					vokhotsk += entry
				else:
					vother += entry
	
	return startstr, endstr, rounded(vokhotsk), rounded(vbering), rounded(vbeaufort), rounded(vchukchi), rounded(vess), rounded(vlaptev), rounded(vkara), rounded(vbarents), rounded(vgreenland), rounded(vcab), rounded(vcaa), rounded(vbaffin), rounded(vhudson), rounded(vother), rounded(vtotal), rounded(volume_uncertainty)#, rounded(area), rounded(extent)	

def createAverage(date):
	startyear = 2014 if date.month <= 4 else 2013
	anomyears = 10
	multiplier = 1.0/anomyears 
	
	folderMonth = date.month
	dayOfYear = date.timetuple().tm_yday
	
	griddedThickness = getGriddedThickness(datetime(startyear, date.month, date.day))
	landmask = insertCryosatDataInNsidcMask(griddedThickness, dayOfYear, startyear, dummyvalue)*multiplier
	print('plotting cryosat anomaly', date)
	
	for k in range(anomyears-1):
		compyear = startyear + 1 + k
		griddedThickness = getGriddedThickness(datetime(compyear, date.month, date.day))
		maskbis = insertCryosatDataInNsidcMask(griddedThickness, dayOfYear, compyear, dummyvalue)
		landmask = addMasks(landmask, maskbis, multiplier, dummyvalue)
	
	savedFileName = 'data/avg/cryosat-smos-avg-' + str(startyear) + '-to-' + str(startyear + 9) + '-' + padzeros(date.month) + padzeros(date.day) + '.csv'
	landmask = np.round(1000*landmask)
	np.savetxt(savedFileName, landmask, delimiter = ',', fmt='%d') #, fmt="%.3f", fmt='%d'

def plotDate(date):
	griddedThickness = getGriddedThickness(date)

	dayOfYear = date.timetuple().tm_yday
	multiplier = 1
	landmask = insertCryosatDataInNsidcMask(griddedThickness, dayOfYear, date.year, dummyvalue)
	landmask = interpolate(landmask, dummyvalue, False)

	plotTitle = "CryoSat-SMOS sea ice thickness " + str(date.day) + " " + monthNames[date.month-1] + " " + str(date.year)
	filename = 'cryosat-smos-thickness-' + str(date.year) + padzeros(date.month) + padzeros(date.day)
	dropboxFilename = ''			
	plotThickness(landmask,plotTitle,filename,dropboxFilename)
			
def uploadToGoogleDrive():
	upload_to_google_drive.replace_file_in_google_drive('1jSihYCk2KkQuMvw1TAldinJy5WGLdygQ','cryosat-smos-volume-cab.png')
	upload_to_google_drive.replace_file_in_google_drive('1477yE9AJBPcH8Pz7QZA21Fjj9g8ipKVg', "cryosat-smos-volume-caa.png")
	upload_to_google_drive.replace_file_in_google_drive('1DON43_2oHN4T8yvpm4xV49ONZmFZ7mvw',"cryosat-smos-volume-beaufort.png")
	upload_to_google_drive.replace_file_in_google_drive('1MrWT5RScXMmojfJFmOKo8P0TXQelGolp',"cryosat-smos-volume-chukchi.png")
	upload_to_google_drive.replace_file_in_google_drive('1_B4Ylz7H4FA5ezO70fOf8HjJn1Fl1FW8',"cryosat-smos-volume-bering.png")
	upload_to_google_drive.replace_file_in_google_drive('13gWz-I_ya0mDsi-OGTInKKuitkBmTvXk',"cryosat-smos-volume-ess.png")
	upload_to_google_drive.replace_file_in_google_drive('1_djGUUYmXCp9nmHwoxQTfLFlSNzHFDLV',"cryosat-smos-volume-laptev.png")
	upload_to_google_drive.replace_file_in_google_drive('1bU5B2oD_IhsCiFqycAMNPgz6h9n-1j8O',"cryosat-smos-volume-kara.png")
	upload_to_google_drive.replace_file_in_google_drive('1WNW-kRoDeoa3Nohejq56ECqZmwqVmxda',"cryosat-smos-volume-barents.png")
	upload_to_google_drive.replace_file_in_google_drive('1GeO_esf9Dw98tD1gK4uxNHTnuEOe4KXW',"cryosat-smos-volume-greenland.png")
	upload_to_google_drive.replace_file_in_google_drive('1b0igRzMVHUqynwEIFmYt6kFWdFLkc2Ji',"cryosat-smos-volume-baffin.png")
	upload_to_google_drive.replace_file_in_google_drive('1fzHE-S8sC_p3IqFkOKQBGZFXGKsObvqW',"cryosat-smos-volume-hudson.png")
	upload_to_google_drive.replace_file_in_google_drive('15jjBCAVOFWzOzDQeTLoyHntXyD328ZVL',"cryosat-smos-volume-okhotsk.png")

mask = np.loadtxt(open("regional-mask.csv", "rb"), delimiter=",", skiprows=0)
refmask = np.loadtxt(open("analysis_sea_ice_thickness_20220415.csv", "rb"), delimiter=",", skiprows=0)
lat = np.loadtxt(open("lat.csv", "rb"), delimiter=",", skiprows=0)
lon = np.loadtxt(open("lon.csv", "rb"), delimiter=",", skiprows=0)

"""
Alternative way to load the regional mask:

maskfilename = "./data/awi-siral-l3c-sithick-cryosat2-rep-nh_25km_ease2-202204-fv2p5[1].nc"
file = Dataset(maskfilename, 'r', format="NETCDF4")
mask = file.variables['region_code'][:].squeeze()
file.close()
"""

dummyvalue=10
thicknessmax = 4.0
anomalymax = 1.0
anomyears = 10 # 10 years in anomaly base

auto = True

datadir = 'LATEST'
filename = ''
if (len(sys.argv) > 1) :
	datadir = sys.argv[1]
if (len(sys.argv) > 2) :
	filename = sys.argv[2]

dummyvalue=10

thicknessmax = 4.0
anomalymax = 1.0
anomyears = 10 # 10 years in anomaly base

if auto:
	plotCryosatThickness = True
	plotCryosatAnomaly = True

	latestDate = downloadNewFiles()

	date = latestDate
	date = date - timedelta(days = 3)
	griddedThickness = getGriddedThickness(date)
		
	dayOfYear = date.timetuple().tm_yday
	print('day',date.day)

if plotCryosatThickness:
	multiplier = 1
	landmask = insertCryosatDataInNsidcMask(griddedThickness, dayOfYear, date.year, dummyvalue)
	landmask = interpolate(landmask, dummyvalue, False)

	plotTitle = "CryoSat-SMOS sea ice thickness " + str(date.day) + " " + monthNames[date.month-1] + " " + str(date.year)
	filename = 'cryosat-smos-thickness-' + str(date.year) + padzeros(date.month) + padzeros(date.day)
	dropboxFilename = 'cryosat-smos-thickness-latest'			
	plotThickness(landmask,plotTitle,filename,dropboxFilename)

if plotCryosatAnomaly:
	multiplier = -1.0/anomyears 
	landmask = insertCryosatDataInNsidcMask(griddedThickness, dayOfYear, date.year, dummyvalue)
	print('plotting cryosat anomaly', date)
	folderMonth = date.month
	for k in range(anomyears):
		compyear = date.year - k - 1 # year-k-1
		griddedThickness = getGriddedThickness(datetime(compyear,date.month,date.day))
		maskbis = insertCryosatDataInNsidcMask(griddedThickness, dayOfYear, compyear, dummyvalue)
		landmask = addMasks(landmask, maskbis, multiplier, dummyvalue)

	landmask = interpolate(landmask, dummyvalue, True)

	print(date.day)
	plotTitle = "CryoSat-SMOS thickness anomaly " + str(date.day) + " " + monthNames[date.month-1] + " " + str(date.year) + " vs " + str(date.year-10) + "-" + str(date.year-1)
	filename = 'cryosat-smos-thickness-anomaly-' + str(date.year) + padzeros(date.month) + padzeros(date.day)
	dropboxFilename = 'cryosat-smos-thickness-anomaly-latest'
	plotAnomaly(landmask,plotTitle,filename,dropboxFilename)

if auto:
	animationFileName = 'animation_cryosat_smos_latest.gif' 
	frames = 10
	for k in range(frames):
		previousdate = date - timedelta(days = k)
		filename = 'cryosat-smos-thickness-' + str(previousdate.year) + padzeros(previousdate.month) + padzeros(previousdate.day) + '.png'
		if not os.path.isfile(filename):
			plotDate(previousdate)
		
	make_animation.makeAnimation(date, frames, animationFileName, lambda date: 'cryosat-smos-thickness-' + str(date.year) + padzeros(date.month) + padzeros(date.day) + '.png')
	if putOnDropbox:
		dropbox_client.uploadToDropbox([animationFileName])

	outFile = open('cryosat-smos-regional-volume.csv', 'a', newline='')
	csvFile = csv.writer(outFile)

	time.sleep(3)
	regional_python_graphs.plotRegionalGraphs()

	time.sleep(3)
	uploadToGoogleDrive()