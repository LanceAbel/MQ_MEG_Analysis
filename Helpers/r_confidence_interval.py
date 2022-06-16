import helpers
from helpers import *
from config import *
from participant_data import *
from participant import Participant

import mne
from mne import io
from mne.stats import permutation_cluster_test
from mne.datasets import sample
import glob
import numpy as np
import pandas as pd
import scipy.io
from scipy.stats.stats import pearsonr
import math
from autoreject import get_rejection_threshold
from autoreject import Ransac
from autoreject import (AutoReject, set_matplotlib_defaults)  # noqa
from autoreject import get_rejection_threshold  # noqa
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mne.preprocessing import ICA
from mne.preprocessing import find_bad_channels_maxwell
from mne.viz import plot_topomap
from mne.viz import plot_compare_evokeds

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


CLUSTER_CUTOFF             = 0.05 # p-value cutoff


# 20% of 1730 trials in each participant, for 72 participants
print(0.2*1730*72)
# Actual number was 24545 when I ran the whole pipeline
n = 24545
z_null = 0
z_critical = st.norm.ppf(1 - CLUSTER_CUTOFF*0.5) # 1.96 at p=0.05
print(z_critical)
z_confidence = [z_null-z_critical*math.sqrt(1/(n-3)),z_null+z_critical*math.sqrt(1/(n-3))]
print(z_confidence)
r_confidence = [fisher_z_to_r(z_confidence[0]), fisher_z_to_r(z_confidence[1])]
print(r_confidence)





#Bonferroni adjustment:
CLUSTER_CUTOFF = CLUSTER_CUTOFF / (4*13)
z_critical = st.norm.ppf(1 - CLUSTER_CUTOFF*0.5) # 1.96 at p=0.05
z_confidence = [z_null-z_critical*math.sqrt(1/(n-3)),z_null+z_critical*math.sqrt(1/(n-3))]
print(z_confidence)
r_confidence = [fisher_z_to_r(z_confidence[0]), fisher_z_to_r(z_confidence[1])]
print(r_confidence)


# Compute z to compare to critical z
n1 = 13099
n2 = 24625
zr1  = -0.145
zr2 = -0.086
z = (zr1-zr2)/math.sqrt((1/(n1-3))+(1/(n2-3)))
print(z)

# Compute z to compare to critical z
n1 = 24545
n2 = 24545
zr1  = -0.033
zr2 = -0.072
z = (zr1-zr2)/math.sqrt((1/(n1-3))+(1/(n2-3)))
print(z)



compare_two_corrs(0.012,0,n1,n2)
compare_two_corrs(0.033,0,n1,n2)