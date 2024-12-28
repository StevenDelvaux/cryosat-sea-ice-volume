from decouple import config
import dropbox

def uploadToDropbox(filenames, folder = ''):
	client = getClient()	
	for computer_path in filenames:
		print("[UPLOADING] {}".format(computer_path))
		dropbox_path= "/" + computer_path
		client.files_upload(open(folder + computer_path, "rb").read(), dropbox_path, mode=dropbox.files.WriteMode.overwrite)
		print("[UPLOADED] {}".format(computer_path))
		
def downloadFromDropbox(filenames):
	client = getClient()
	for dropbox_path in filenames:
		print("[DWONLOADING] {}".format(dropbox_path))
		computer_path= dropbox_path
		with open(computer_path, 'wb') as f:
			metadata, res = client.files_download(path= "/" + dropbox_path)
			f.write(res.content)
		print("[DOWNLOADED] {}".format(dropbox_path))
		
def getClient():
	dropbox_access_token = config('DROPBOX_ACCESS_TOKEN')
	app_key = config('APP_KEY')
	app_secret = config('APP_SECRET')
	oauth2_refresh_token = config('OAUTH2_REFRESH_TOKEN')
	client = dropbox.Dropbox(oauth2_access_token=dropbox_access_token,app_key=app_key,app_secret=app_secret,oauth2_refresh_token=oauth2_refresh_token)
	print("[SUCCESS] dropbox account linked")
	return client