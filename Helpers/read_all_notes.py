# Read all notes taken by previous researcher

import time
from time import sleep
import glob

# Has to be run in jupyter notebook, the sleep() doesn't work
def visitfile(file):
    global all_files
    all_files.append(file)

base = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\\'
import os
file_names = [os.path.join(dp, f) for dp, dn, filenames in os.walk(base) for f in filenames if os.path.splitext(f)[1] == '.txt']

notes_files = []
for file in file_names:
	#print(file)
	if ("NOTES" in file or ".rtf" in file) and "PIPE" not in file:
		print(file)
		with open(file) as f:
			lines = f.readlines()
		print(lines)
		sleep(10)
