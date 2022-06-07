
## Check if the participants' sex and left/right handedness information is present


import os
import sys
from os import listdir
from os.path import isfile, join
import numpy as np
import math
import send2trash
from stat import *

adult = True  # Run on adult subfolders if True, else Child subfolders

# SUBJECTS = ['2748'] # If > 9000, the same predictor sequence would be generated as the tones played will be the same

# Where to look for the sequence of tones
CHILD_BASE = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Child_MEG\\'
ADULT_BASE = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\\'
BASE = ADULT_BASE if adult else CHILD_BASE


def get_data(LOGFILE):
    '''Get data out of file as lines'''
    #try:
    with open(LOGFILE) as myfile:
        content = myfile.readlines()
    return content

def walktree(top, callback):
    '''recursively descend the directory tree rooted at top,
       calling the callback function for each regular file'''

    for f in os.listdir(top):
        pathname = os.path.join(top, f)
        mode = os.stat(pathname)[ST_MODE]
        if S_ISDIR(mode):
            # It's a directory, recurse into it
            walktree(pathname, callback)
        elif S_ISREG(mode):
            # It's a file, call the callback function
            callback(pathname)
        else:
            # Unknown file type, print a message
            #print('Skipping %s' % pathname)
            pass 

def visitfile(file):
    global all_files
    all_files.append(file)

def find_all_files(folder):
    global all_files
    all_files = []
    walktree(folder, visitfile)
    return all_files


# Adults
PROBLEM_PARTICIPANTS   =  ['2607']              # Error loading RAW, also there are up to 11 repetitions here, not 7
#PROBLEM_PARTICIPANTS   += ['2589','2480']      # .CON FILE ISN'T ~15 MINS (suspicious of what was shown, can check to add these back)  
#PROBLEM_PARTICIPANTS   += ['2689']             # No sound detected on Matlab pipeline. But what about in this python pipeline?
# Children
PROBLEM_PARTICIPANTS   += ['2724','2810']       # No events on any sound channel that I can see.


os.chdir(CHILD_BASE)
child_participant_strings = []
base_contents = os.listdir('.')
for item in base_contents:
    if os.path.isdir(item):
        if len(item) == 4: # All participants have a four-number code
            try:
                participant_number = int(item)
                child_participant_strings.append(item)
            except:
                print("Not a participant: ", item )
child_participant_strings = sorted(list(set(child_participant_strings)-set(PROBLEM_PARTICIPANTS)))
print("Children: ", child_participant_strings)


# Set up adult participants
os.chdir(ADULT_BASE)
adult_participant_strings = []
base_contents = os.listdir('.')
for item in base_contents:
    if os.path.isdir(item):
        if len(item) == 4: # All participants have a four-number code
            try:
                participant_number = int(item)
                adult_participant_strings.append(item)
            except:
                print("Not a participant: ", item )

adult_participant_strings = sorted(list(set(adult_participant_strings)-set(PROBLEM_PARTICIPANTS)))
print("Adults: ", adult_participant_strings)






strings = adult_participant_strings if adult else child_participant_strings
failed_to_run = []
if 'SUBJECTS' in locals():
    strings = SUBJECTS

sex_missing = []
handedness_missing = []
print(strings)

for ptcp in strings:
    participant_path = BASE+ptcp
    os.chdir(participant_path)
    try:
        sex = get_data('Sex.txt')
    except:
        sex_missing.append(ptcp)    
    try:
        handedness = get_data('Handedness.txt')
    except:
        handedness_missing.append(ptcp)    
    print(ptcp,sex)

print("Missing sex info for : ", sex_missing)
print("Missing handedness info for : ", handedness_missing)    

# onlyfolders = [f for f in listdir(participant_path) if not isfile(join(participant_path, f))]
# print(onlyfolders)