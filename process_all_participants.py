# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import sys
import io
# Need this when running in sublime
sys.stdout = io.TextIOWrapper(sys.stdout.detach(), encoding = 'utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.detach(), encoding = 'utf-8')



import mne
from mne import io
from mne.stats import permutation_cluster_test
from mne.datasets import sample
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

is_pc = os.name=='nt' == True
SCRIPT_ROOT  = os.getcwd() # Directory this is running from
HELPER_DIR   = SCRIPT_ROOT+"\\Helpers\\"
START_NOTEBOOK_TIME = time.time()

# # Import Helpers
import helpers
# Remove manually-defined bad channels (per participant)
# Remove participants that were excluded from Hannah's experiment (ID < 9000) and mine (>9000)
# %cd E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Lance\MEG\Python MEG
# os.chdir(THIS_DIR) 
from config import *
print("Events to tag: ", events_to_tag)
print("Events to tag rerun: ", events_to_tag_rerun)
print("Conds to compare ", conds_to_compare[events_to_tag[1]])
print("Conds to compare ", conds_to_compare[events_to_tag_rerun[1]])
# ## Set up participants
from participant_data import *
# # PARTICIPANT
os.chdir(SCRIPT_ROOT) 
from participant import *
# # EXPERIMENT
from experiment import *
print(attrs_to_delete)



# ### Add many participants from scratch
KEEP_IN_MEMORY    = False # If true, will keep in memory each Participant object
SAVE_PARTICIPANTS = True

children_memory   = []
adults_memory     = []

MMR_Experiment = Experiment(experiment_numbers=experiment_numbers["Children"], condition = conds_to_compare[events_to_tag[1]][1])
print("Condition of high surprise: ", conds_to_compare[events_to_tag[1]][1])
#child_participant_strings = child_participant_strings[0:5]
#adult_participant_strings = adult_participant_strings[0:5]
print("%s child system, strings: %s"%(str(len(child_participant_strings)), str(child_participant_strings)))
print("%s adult system, strings: %s"%(str(len(adult_participant_strings)), str(adult_participant_strings)))



# # Add participants to the experiment by running pipeline
OLDFILE_CUTOFF_DAYS = 0.001 # Only run the pipeline if the saved file was created > this many days ago. Useful for re-running analysis on all participants after a change
file_num_start  = 0
file_num_end    = 12
child_participant_strings = ['2629'] # sorted(child_participant_strings[min(len(child_participant_strings),file_num_start):min(len(child_participant_strings),file_num_end)],reverse=True)
adult_participant_strings = ['2552'] # sorted(adult_participant_strings[min(len(adult_participant_strings),file_num_start):min(len(adult_participant_strings),file_num_end)],reverse=True)

# Add participants to the experiment by running pipeline
children_additions = {"successful": [], "errors": []}
adult_additions    = {"successful": [], "errors": []}
evoked_trigger_leakages = []

timestamp_file = r'C:\Users\Lance\Desktop\CanDelete\\Temp.txt'
with open(timestamp_file, 'w') as f:
    f.write('Temp file')
time_now = os.path.getmtime(timestamp_file)


exception_log = []
all_warnings = []
r = 1

    

for adult_string in adult_participant_strings:
    print("@@@@@@@@@@@@@ Participant # ", r, " of", len(child_participant_strings)+len(adult_participant_strings), adult_string)   
    num_calls = {}    
    try: 
        file_names = glob.glob(f"{SAVE_DIR}*.pickle")
        file_found = False
        days_since = 9999        
        for file in file_names:
            if adult_string in file:
                file_found = True
                time_mod = os.path.getmtime(file)
                days_since = (time_now - time_mod)/60/60/24

        if not file_found or days_since > OLDFILE_CUTOFF_DAYS: # File wasn't found or is old
            print(adult_string, days_since, " Days")
            ptcp = Participant(is_adult_folder=True, p_id=adult_string, experiment_number=1)
            ptcp.basic_cleaning()
            ptcp.initiate_events()
            ptcp.downsample()
            ptcp.more_processing()
            MMR_Experiment.update_tracking(ptcp)
            #process_trigger_leakage(ptcp)
            #evoked_trigger_leakages.append(ptcp.evoked_trigger_leak)
            del_attributes(ptcp, attrs_to_delete)
            if SAVE_PARTICIPANTS:
                ptcp.save_participant()
            adult_additions['successful'] = adult_additions['successful'] + [adult_string]
            if KEEP_IN_MEMORY: #  or r == 1
                adults_memory.append(ptcp)     
            gc.collect()
            clear_output()
            print("Warnings were: ", ptcp.warnings)
            all_warnings.append(ptcp.warnings)
    except Exception as e:
        print("############## PROBLEMS ADDING ADULT ###############")
        print(adult_string, e)
        exception_log.append([adult_string, e])        
        adult_additions['errors'] = adult_additions['errors'] + [adult_string]
    print("Exception log: ", exception_log)
    r+=1

for child_string in child_participant_strings:
    print("@@@@@@@@@@@@@ Participant # ", r, " of", len(child_participant_strings)+len(adult_participant_strings), child_string)    
    num_calls = {}       
    try:
        file_names = glob.glob(f"{SAVE_DIR}*.pickle")
        file_found = False
        days_since = 9999
        for file in file_names:
            if child_string in file:
                file_found = True
                time_mod = os.path.getmtime(file)
                days_since = (time_now - time_mod)/60/60/24

        if not file_found or days_since > OLDFILE_CUTOFF_DAYS:
            print(child_string, days_since, " Days")
            #if not file_found:        
            ptcp = Participant(is_adult_folder=False, p_id=child_string, experiment_number=1)
            ptcp.basic_cleaning()
            ptcp.initiate_events()
            ptcp.downsample()
            ptcp.more_processing()
            MMR_Experiment.update_tracking(ptcp)
            #process_trigger_leakage(ptcp)
            #evoked_trigger_leakages.append(ptcp.evoked_trigger_leak)
            del_attributes(ptcp, attrs_to_delete)
            if SAVE_PARTICIPANTS:
                ptcp.save_participant()
            children_additions['successful'] = children_additions['successful'] + [child_string]
            if KEEP_IN_MEMORY: # or r == 1
                children_memory.append(ptcp)
            gc.collect()
            clear_output()
            print("Warnings were: ", ptcp.warnings)
            all_warnings.append(ptcp.warnings)
    except Exception as e:
        print("############## PROBLEMS ADDING CHILD ###############")
        print(child_string, e)
        exception_log.append([child_string,e])
        children_additions['errors'] = children_additions['errors'] + [child_string]
    print("Exception log: ", exception_log)
    r+=1
 

print(" @@@@@@@@@@@@ FINISHED ADDING @@@@@@@@@@@@")    
print(children_additions)
print(adult_additions)

print("\n")
print("Exceptions: ", exception_log)
print(time.time(), "Took ", (time.time()-START_NOTEBOOK_TIME)/60, " Minutes")


print("Printing all warnings:")
if len(all_warnings) > 0:
    for ptcp_warns in all_warnings:
        print(ptcp_warns)