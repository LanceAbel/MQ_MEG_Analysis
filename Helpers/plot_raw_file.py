


import sys
import io
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

is_pc = os.name=='nt' == True



### FROM SINGLE PICKLE FILE
class Participant():
    def __init__(self, is_adult_folder, p_id=None, experiment_number=1):
        
        pass


## Compare adult vs child system recording
raw_adult = mne.io.read_raw_kit(
    input_fname=r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\9000\\9000_JW_21_11_19_experiment1.con',
    slope="+",
    stim_code="channel",
    stimthresh=2,
    preload=True,
    allow_unknown_format=False,
    verbose=False)

raw_child = mne.io.read_raw_kit(
    input_fname=r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Child_MEG\9001\\9001_AH_21_12_15_experiment1.con',
    slope="+",
    stim_code="channel",
    stimthresh=2,
    preload=True,
    allow_unknown_format=False,
    verbose=False)


DURATION 	= 10 # (secs)
FREQ 		= 1000  # /sec
N_CHANNELS 	= 1
## See that the plots line up
START_ADULT = 122
DELAY_CHILD = 3.0 # How long after the adult system started recording I started on Child (secs)
START_CHILD = START_ADULT - DELAY_CHILD
raw_adult.plot(start=START_ADULT, duration=DURATION,n_channels=N_CHANNELS,remove_dc=False)
plt.show()

raw_child.plot(start=START_ADULT-DELAY_CHILD, duration=DURATION,n_channels=N_CHANNELS,remove_dc=False)
plt.show()







