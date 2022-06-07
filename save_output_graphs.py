# Runs all statistical analyses for each condition/model used, saves outputs (graphs of GFPs, ERFs, significant temporal clusters, significant spatiotemporal clusters, regression analysis outputs etc)
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
start_time = time.time()
from IPython.core.display import display, HTML
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



fontsize = 11 # Font size for titles, where applicable
percentiles = 10, 20
age_bounds_low=[3,5]
age_bounds_high=[10,20]
EPOCH_START_ANALYSIS, EPOCH_END_ANALYSIS = 0.0, 0.295 # 0.295 #  # 0.13, 0.185 
CONSTANT_RADIUS_MULT = 1 # 0.85
EXP_BASE          = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Experiments\\'
EXP_BASE          = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\ExperimentsNew\\'
os.chdir(EXP_BASE)
base_contents = os.listdir('.')
STAT_TO_USE = 'corr' # 'corr'
# Extracts the evoked for each condition tested and the group average
CHILD_BASE          = BASE_FOLDER+'\Child_MEG\\'
ADULT_BASE          = BASE_FOLDER+'\Adult_MEG\\'
child_participant_strings, adult_participant_strings = set_up_participants() # > set_up_participants found in experiment.py
print("Number children", len(child_participant_strings), "Number adults ", len(adult_participant_strings))
print(child_participant_strings)
print(adult_participant_strings)



import scipy.stats as st
base_contents_to_use = base_contents[-1:]
start_time = time.time()
problems = []
MMR_Experiment = Experiment(experiment_numbers=experiment_numbers["Children"], condition = conds_to_compare[events_to_tag[1]][1])
MMR_Experiment.setup_to_run()
for folder in base_contents_to_use:
    base = EXP_BASE+folder
    if os.path.isdir(base):
        tgt_dir = base +"\\" # +'\\Standardised\\'    
        print(tgt_dir)

        #clear_output()
        MMR_Experiment.all_stats_loop(tgt_dir, MMR_Experiment=None)

        # Also produce files post-standardisation
        tgt_dir = base +"\\" +'\\Standardised\\'    
        print(tgt_dir)        
        MMR_Experiment.all_stats_loop(tgt_dir, MMR_Experiment='NotNone')

    clear_output()
print("Problems with", problems)
print("Took", (time.time()-start_time)/60, " minutes ")