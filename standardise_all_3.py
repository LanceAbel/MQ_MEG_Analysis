# Standardise the data contained in all the participant files within a folder
import mne
from mne import io
from mne.stats import permutation_cluster_test

import glob
import numpy as np
import pandas as pd
import scipy.io
import math

from autoreject import get_rejection_threshold
from autoreject import Ransac
from autoreject import (AutoReject, set_matplotlib_defaults)  # noqa
from autoreject import get_rejection_threshold  # noqa
from matplotlib import pyplot as plt
#import Image
import matplotlib.pyplot as plt
from mne.preprocessing import ICA
from mne.preprocessing import find_bad_channels_maxwell
import gc
import os
import time
import re

import pickle
import hickle
import copy

from pathlib import Path
from IPython.display import clear_output
import inspect
os.chdir(r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Lance\MEG\Python MEG\\')
import numpy as np
import glob
from helpers import *
import pandas as pd
from participant_data import *
from config import *
os.chdir(r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Lance\MEG\Python MEG\\')
from experiment import *
from participant import *




# Extracts the evoked for each condition tested and the group average
CHILD_BASE          = BASE_FOLDER+'\Child_MEG\\'
ADULT_BASE          = BASE_FOLDER+'\Adult_MEG\\'


def apply_z_score_new(ptcp,STANDARDISE=True,channel_wise=False):
	if STANDARDISE:
		for attrib in ['epochs','epochs_ransac','epochs_ransac_autoreject','epochs_ransac_autoreject_unstandardised']:
			if hasattr(ptcp,attrib):
				obj = getattr(ptcp, attrib)
				if channel_wise:
					if type(obj) == type({}):
						keys = list(obj.keys())
						for key in keys:
							print("Applying, ", key)
							obj[key].apply_function(fun=calc_z_score, n_jobs=4)
					else:
						obj.apply_function(fun=calc_z_score, n_jobs=4)
				else:
					try:
						obj.apply_function(fun=calc_z_score)  
						print("H1")                    
					except Exception as e:
						print("H2")
		ptcp.signals_of_interest()	# -> This creates 'evoked_generic','evoked_all']
		ptcp.applied_z_score = True
	else:
		ptcp.applied_z_score = False







EXP_BASE          = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Experiments\\'
EXP_BASE          = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\ExperimentsNew\\' # PE2 Vanilla 20pct -100 to 300\\'
EXP_BASE          = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Experiments Final\\'

os.chdir(EXP_BASE)
print(EXP_BASE)
base_contents = os.listdir('.')
print(base_contents)

problem_folders = []
error_messages = []
folder_num = 1
#for folder in base_contents[6:7]+base_contents[4:5]+base_contents[22:23]:
for folder in base_contents[2:3]: # + base_contents[2:3]: # +  # base_contents[14:15]+base_contents[20:21]: # 
	r = 0
	tgt_dir = EXP_BASE+folder+'\\'
	print("Target dir, ", tgt_dir, " folder #: ", folder_num)
	os.chdir(tgt_dir)
	file_names = glob.glob(f"{tgt_dir}*.pickle")


	#try:
	for file in file_names:
		include = False
		for child_string in child_participant_strings:
			if child_string in file:
				include = True
		for adult_string in adult_participant_strings:
			if adult_string in file:
				include = True            
		if include:

			# Figure out which experiment it is for
			with open(file, "rb") as f:
				ptcp = pickle.load(f)
			print(" folder #: ", folder_num, r, ptcp.p_id, "Target dir, ", tgt_dir)
			condition_to_compare = ptcp.cond_B
			events_to_tag = ptcp.events_to_tag
			#print(ptcp.evoked_all.to_data_frame().head(2))
			#apply_z_score(ptcp)
			apply_z_score_new(ptcp)
			#print(ptcp.evoked_all.to_data_frame().head(2))


			if ("unstandardised" in file or "STD=F" in file) and (" standardised" not in file and "STD=T" not in file): # They're all unstandardised files
				if "STD=F" in file:
					new_file_name = file.replace("STD=F", "STD=T")	
				else:
					new_file_name = file.replace("unstandardised.pickle", "STD=T.pickle")

				pickle.dump(ptcp, open(new_file_name, "wb")) 
			else:
				if (" standardised" not in file and "STD=T" not in file): # And we don't already have a file
					if "STD=F" in file:
						new_file_name = file.replace("STD=F", "STD=T")		
					else:
						new_file_name = file.replace(".pickle", " STD=T.pickle")					
					print(new_file_name)
					pickle.dump(ptcp, open(new_file_name, "wb")) 
		r+=1

	# except Exception as e:
	# 	print("Error ", folder_num, e)


	folder_num+=1
