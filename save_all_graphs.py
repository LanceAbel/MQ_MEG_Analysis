#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Runs all statistical analyses for each condition/model used, saves outputs (graphs of GFPs, ERFs, significant temporal clusters, significant spatiotemporal clusters, regression analysis outputs etc)


# In[2]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
__author__ = "Lance Abel"
__copyright__ = "2022"
__credits__ = "With input from Paul Sowman"
__license__ = "GPL"
__version__ = "27"
__maintainer__ = "Lance Abel"
__email__ = "lance.abel@hdr.mq.edu.au"
__status__ = 

# mne 0.23.4
# numpy 1.19.5
# pandas # 0.25.3
# matplot 3.3.4
# pimpler 1.0.1
# regex 2022.1.18
# scipy 1.5.4
# hickle 4.0.4
# regex 2022.3.15
# autoreject 0.2.2
# pathlib 1.0.1
# joblib 1.1.0
# sklearn 0.24.2
# psutil 5.9.0

"""
#%% required packages
# ! pip install -- upgrade mne
# ! pip install https://api.github.com/repos/autoreject/autoreject/zipball/master
# ! pip install fooof
# ! pip install --upgrade numpy
# ! pip install pyriemann
# !pip install eeglabio
# pip install dss


# In[3]:


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
display(HTML("<style>.container { width:100% !important; }</style>"))
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.options.display.float_format = '{:,.4f}'.format


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


EPOCH_START_ANALYSIS, EPOCH_END_ANALYSIS = 0.005, 0.295 # 0.295 #  # 0.13, 0.185 

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

base_contents
base_contents_to_use = base_contents[-9:-8] #not done
base_contents_to_use

start_time = time.time()
problems = []
MMR_Experiment = Experiment(experiment_numbers=experiment_numbers["Children"], condition = conds_to_compare[events_to_tag[1]][1])
MMR_Experiment.setup_to_run()
for folder in base_contents_to_use:
    base = EXP_BASE+folder
    if os.path.isdir(base):
        tgt_dir = base +"\\" # +'\\Standardised\\'    
        print(tgt_dir)

#         #clear_output()
#         MMR_Experiment.all_stats_loop(tgt_dir, age_group="Both", MMR_Experiment=None)

        #Also produce files post-standardisation
        tgt_dir = base +"\\" +'\\Standardised\\'    
        print(tgt_dir)        
        MMR_Experiment.all_stats_loop(tgt_dir, age_group="Both", MMR_Experiment='NotNone')

    clear_output()
print("Problems with", problems)



len(MMR_Experiment.participants)


print("Took", (time.time()-start_time)/60, " minutes ")


