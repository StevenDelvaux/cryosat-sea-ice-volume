from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from googleapiclient.discovery import build
import os
import os.path
from googleapiclient.http import MediaFileUpload
from decouple import config

def get_credentials(SCOPES):
	creds = None
	credentials_filename = "token.json"
	if not os.path.exists(credentials_filename):
		google_drive_credentials = config('GOOGLE_DRIVE_CREDENTIALS')
		with open(credentials_filename, "w") as local_file:
			local_file.write(google_drive_credentials)
	
	creds = Credentials.from_authorized_user_file(credentials_filename, SCOPES)
	
	if not creds or not creds.valid:
		print('credentials invalid')
		if creds and creds.expired and creds.refresh_token:
			print('credentials try token refresh')
			creds.refresh(Request())
			with open(credentials_filename, 'w') as token:
				token.write(creds.to_json())
	return creds

def replace_file_in_google_drive(file_id,local_path):
	credentials = get_credentials(["https://www.googleapis.com/auth/drive.file"])
	drive_service = build('drive', 'v3', credentials=credentials)

	file_metadata = {'name': local_path}
	media = MediaFileUpload(local_path, mimetype='image/png')
	file = drive_service.files().update(fileId = file_id, media_body=media).execute()
	print(F'File ID: {file.get("id")}')
