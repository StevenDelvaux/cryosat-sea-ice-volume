def getLastSavedDay(filename):
	print('inside last saved day', filename)
	with open(filename, 'r') as f:
		lastline = f.readlines()[-1]
	splitted = lastline.split(",", 3)
	print('inside last saved day', ' last line ', lastline)
	return int(splitted[0]), int(splitted[1])
