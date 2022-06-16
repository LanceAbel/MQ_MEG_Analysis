import mne
from mne import io
from mne.stats import permutation_cluster_test
from mne.datasets import sample
import glob
import numpy as np
import pandas as pd
import scipy.io
import math
from matplotlib import pyplot as plt
from autoreject import get_rejection_threshold
from autoreject import Ransac
from autoreject import (AutoReject, set_matplotlib_defaults)  # noqa
from autoreject import get_rejection_threshold  # noqa
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
import scipy.stats as st

BASE_FOLDER			= r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\\' # Data root folder before \Adult & \Child
SAVE_DIR			= BASE_FOLDER+'\ExperimentsNew\\' # Directory to save (processed) experimental data in

from participant_data import *
from helpers import *


LOGGING				       = True       # Each function will track key points, to identify where Exceptions occur
STANDARDISE                = False      # Whether to standardise the ERF values or not

#### CONFIG 1 - Stats parameters
#Time window of analysis
EPOCH_START                = -0.1
EPOCH_END                  = 0.3
EPOCH_START_ANALYSIS       = 0.0       # The min latency at which we expect a MMR based on previous research (secs)
EPOCH_END_ANALYSIS         = 0.295       # The max latency at which we expect a MMR based on previous research (secs)
EPOCH_SIZE_MS              = 5
NUM_BINS_SURPRISE          = 5          # x binsA, we divide the sample into the top and bottom 100%/x bins as high and low suirprise respectively. So 5 gives 0-20th and 80-100th percentiles of surprise
SEC_TO_MS                  = 1000
EPSILON_TIME               = 0.01       # Extra bit added so that quantities rounded down map to the closest integer

NUM_EPOCH_SAMPLES          = int((EPOCH_END_ANALYSIS-EPOCH_START_ANALYSIS)*1000 / EPOCH_SIZE_MS + 1) # Number of time samples in an epoch
CLUSTER_CUTOFF             = 0.05       # p-value cutoff for all significance testing
N_PERMUTATIONS             = 10000      # For permutation testing, how many permutations to test
STAT_TO_USE                = 'corr'     # Which statistic to use when analysing relationship between the predictor and the GFP or ERF



#### CONFIG 2 - preprocessing options
# Data cleaning
RUN_RANSAC                 = False       # Set to False if rerunning a condition i.e. when using a saved file in which Ransac was already run (and in which you now just want to change the condition...)
RUN_AUTOREJECT             = False
#ICA
RUN_ICA                    = False
SAVE_ICA                   = True if RUN_ICA else False
ICA_NUM_COMPONENTS         = 15
MIN_CHANNELS_ICA           = 3          # Min number of channels a component must register as significant in order to remove
ICA_THRESHOLD_EOG          = 0.8
ICA_THRESHOLD_ECG          = 0.95

#Save
SAVE_GROUP_ANALYSIS        = False      # Will save to .hickle/.pickle if True

# CPUs to utilise. Going above 1 occasionally causes some errors on Ransac. Optimum would be try to with 6/8, then if it fails re-do
NUM_CPUS_Ransac            = 1                    
NUM_CPUS_LineFilter        = 1          # 4
NUM_CPUS_Autoreject        = 1          # 6      # Doesn't appear to fail with > 6.
NUM_CPUS_Other             = 1
NUM_CPUS_Stats             = 4

# Filters
HP_FREQ                    = 0.5        # Min freq to keep
LP_FREQ                    = 60         # Max freq to keep
MACHINE_DOWNSAMPLE_FREQ    = 200        # What the MEG machine itself downsampled to
DOWNSAMPLE_FREQ            = 200        # What to downsample to
LINE_FREQ                  = 50


 
AGE_CUTOFF					= 10        # When splitting group into age < X  or age  > X , what to use for X
CUTOFF_FILE_NAMES			= 9000      # I (in 2022) started labelling new participants with IDs 9000 and above

# Note: 2629 , 2687, 2695, 2632 have the opposite of MMN and a lot of HF noise. Certainly not excluding them because of this, just interesting to note.
# Normal ones have lower GFP for deviant, the above have higher GFP for deviant


CONT_PREDICTOR_USE_DIFF_CUTOFF  = True  # If True, will split the predicted surprise not into   >=x and <x ,  but into >=x and <=y where y!=x.
                                        # This choice only matters for predictors that take on non-binary values
                                        # To aid comparability between the HGF and deviants, set this to True

##### Choose predictor, definition of events
# Which predictors we are assessing - pick one
events_to_tag = {0: "frequencies",             # Choose a mode for which events we care about within participant data by setting its value to 1. Can change this later and continue running
                 
                 0: "deviants", 1: "deviants_custom", 0: 'deviants_custom_combine', 0: "deviants_specific",
                 
                 0: "hgf_pe2",   0: "hgf_pe2_mod_baked",    0: "hgf_pe2_mod",    
                 0: "hgf_pwpe2", 0: "hgf_pwpe2_mod_baked",  # 0: "hgf_pwpe2_mod",   # Wasn't created
                 0: "hgf_pe3",   0: "hgf_pe3_mod_baked",    # 0: "hgf_pe3_mod",     # Wasn't created    
                 0: "hgf_pwpe3", 0: "hgf_pwpe3_mod_baked",  # 0: "hgf_pwpe3_mod",   # Wasn't created 


                 0: "PS", 0: "CS", 0: "BS"
                 }

# New conditions to apply (if CHANGE_CONDITION==True)
events_to_tag_rerun = {0: "frequencies",             
                 
                 0: "deviants", 1: "deviants_custom", 0: 'deviants_custom_combine', 0: "deviants_specific",
                 
                 0: "hgf_pe2",   0: "hgf_pe2_mod_baked",        0: "hgf_pe2_mod",
                 0: "hgf_pwpe2", 0: "hgf_pwpe2_mod_baked",      # 0: "hgf_pwpe2_mod",   # Wasn't created
                 0: "hgf_pe3",   0: "hgf_pe3_mod_baked",        # 0: "hgf_pe3_mod",     # Wasn't created    
                 0: "hgf_pwpe3", 0: "hgf_pwpe3_mod_baked",      # 0: "hgf_pwpe3_mod",   # Wasn't created 

                 0: "PS", 0: "CS", 0: "BS"
                 }

# Main comparison conditions - these are just unique codes to identify the (supervised feature) conditions compared
conds_to_compare = {'frequencies': ['600', '700'], # Can change to e.g. 650 vs 800 to compare these
                    
                    'deviants': ['88', '99'], 'deviants_custom': ['887','997'], 'deviants_custom_combine': ['991','999'], # 999 is 5 or more, 991 is 4 or less
                    'deviants_specific': ['994','995'], # e.g. a deviant after 4 repetitions should be much more surprising than one after 5
                    #'deviants_specific': ['997','991'], # e.g. a deviant after 4 repetitions should be much more surprising than one after 5
                   
                    'hgf_pe2': ['11','12'],   'hgf_pe2_mod_baked': ['15','16'],                             'hgf_pe2_mod': ['13','14'], 
                    'hgf_pwpe2': ['21','22'], 'hgf_pwpe2_mod_baked': ['23','24'],                           # 'hgf_pwpe2_mod':  ['13','14'], # Wasn't created 
                    'hgf_pe3': ['31','32'],   'hgf_pe3_mod_baked': ['33','34'],                             # 'hgf_pe3_mod':    ['15','16'],   # Wasn't created 
                    'hgf_pwpe3': ['41','42'], 'hgf_pwpe3_mod_baked': ['43','44'],                           # 'hgf_pwpe3_mod':  ['17','18'], # Wasn't created                  
                    
 
                    'PS': ['71','72'],
                    'CS': ['81','82'],
                    'BS': ['91','92']
                   } # The labels of the event type that we care about. Can change this later and continue running


HGF_TIMES                       = [str(346*x) for x in range(1,6)] + ['2000']  # The length of the sequences used to train various HGFs as predictors 
RELEVANT_LENGTH_PREDICTOR       = HGF_TIMES[-1] # Time to use for main HGF tests...


experiment_numbers = {"Adults": [1], "Children": [1]} # Not implemented yet, but will pick the experiment of interest when multiple recordings have been made from the same individuals








#### CONFIG 3 - less likely to change - stimulus & data characteristics
NUM_TONES                  = 7                   # Number of distinct auditory tones played
MAX_NUM_REPEATS            = 7                   # Number of timese each auditory tone repeats before a 'deviant'
ALL_TONE_FREQUENCIES       = [500+x*50 for x in range(0,NUM_TONES)]
TONE_DURATION              = 0.07                # secs

# Dictionary of exceptions to the rule of which one the audio channel is. Haven't been able to auto-detect this.
SOUND_EXCEPTIONS           = {'3448': 183  # } 
                                # Was 132 for these:
                                # '3374': 135, '3419': 135, '3421': 135, '3422': 135, '3423': 135, '3429': 135,
                                # '3434': 135, '3438': 135, '3439': 135
                              }
for lance_ptcp in [9000, 9001]: # + list(range(9002,9008)
    SOUND_EXCEPTIONS[str(lance_ptcp)] = 183
    
    
SOUND_DELAY_1              = 0.044               # Backup sound delay for participants with ID < SOUND_DELAY_2_BEGINS
SOUND_DELAY_2              = 0.244               # Backup sound delay for participants with ID >= SOUND_DELAY_2_BEGINS
SOUND_DELAY_2_BEGINS       = 9000                # First participant ID to have the longer sound latency(SOUND_DELAY_2)

TRIGGER_EXCEPTIONS         = {}                  # Dictionary of exceptions to the rule of which one the MEG trigger (events) channel is. Haven't been able to auto-detect this.





is_adult = False
# Only needed when running single participants
stim_begin_adult = 194
stim_begin_child = 146
# Needing when running experiment or group (adding in many participants)
exp_num                     = 1
num_meg_chs_adult           = 160                   # How many magnetometer channels will be found in MNE Raw if the system was KIT160 adult system
num_meg_chs_child           = 125                   # How many magnetometer channels will be found in MNE Raw if the system was KIT160 child system
num_chs_adult               = 257                   # How many will be found in MNE Raw if the system was KIT160 adult system
num_chs_child               = 193                   # How many will be found in MNE Raw if the system was KIT160 child system
print("Adult :", is_adult, "Experiment number: ", exp_num, " Number of channels on system adult: ", num_chs_adult-1, " Number of channels child: ", num_chs_child-1)
audio_channel_adult         = 167                   # 167 (MISC 007)
audio_channel_child         = 135                   # 135 (MISC 010)
NUM_CHANS_TO_KEEP           = 124



# Plot settings     
DURATION        = 2
N_CHANS         = 3
START_TIME      = 120
RESOLUTION      = 64 * 3
SIZE            = 2
YLIM            = 20
DPI             = 600     # Plot quality
FONTSIZE_ALL    = 14
FONTSIZE_TITLE  = 19
FONTSIZE_LABELS = 18  # Title on x and y axis
FONTSIZE_AXES   = 16  # Values on x and y axes
FONTSIZE_LEGEND = 13  # Legend
PADDING         = 20  # Space between axis values and axis label


sphere_x, sphere_y, sphere_z, radius = 0.007, -0.017, 0.00, 0.107
CONSTANT_RADIUS_MULT = 1 # 0.85








string_save = events_to_tag[1] + " "
if RUN_RANSAC:
    string_save+="RS=T,"
else:
    string_save+="RS=F,"
if RUN_AUTOREJECT:
    string_save+="AR=T,"
else:
    string_save+="AR=F,"
if RUN_ICA:
    string_save+="ICA=T,"
else:
    string_save+="ICA=F,"
if STANDARDISE:
    string_save+="STD=T,"
else:
    string_save+="STD=F,"
string_save += " " +str(int(1000*EPOCH_START_ANALYSIS)) + "to" +str(int(1000*EPOCH_END_ANALYSIS))
if CONT_PREDICTOR_USE_DIFF_CUTOFF:
    string_save += " "+str(int(100/NUM_BINS_SURPRISE)) + "%"+ "vs"+str(int(100/NUM_BINS_SURPRISE)) + "%"
else:
    string_save +=" " +str(int(100/NUM_BINS_SURPRISE)) + "%"+ "vs"+str(100-int(100/NUM_BINS_SURPRISE)) + "%"
print("File to be saved as :", string_save)
SAVE_EXPERIMENT_NAME = string_save