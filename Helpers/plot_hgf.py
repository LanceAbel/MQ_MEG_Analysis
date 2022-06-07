import sys
import io

import glob
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from mne.preprocessing import ICA
from mne.preprocessing import find_bad_channels_maxwell


base_folder = r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\HGF Settings999\\'
base_folder2 = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Common\\'

def get_data(file_path):
    '''Get data out of file as lines'''
    with open(file_path) as myfile:
        content = myfile.readlines()
    return content
def process_file(lines):
    return [float(x.replace("\n","")) for x in lines]
def process_file2(string):
    string = string[0]
    string = string.split(",")
    string = string[0:-1]
    string = [float(x.replace("\n","")) for x in string]
    return string


file_dict = {   #"HGF PE2": base_folder+'HGF1730 PE2_new.txt',
                #"HGF PE3": base_folder+'HGF1730 PE3_new.txt',
                #"HGF PWPE2": base_folder+'HGF1730 PWPE2_new.txt',
                #"HGF PWPE3": base_folder+'HGF1730 PWPE3_new.txt',
                "HGF Integrated PE2": base_folder+'HGF1730 PE2_mod_baked_new.txt',
                "HGF Integrated PE3": base_folder+'HGF1730 PE3_mod_baked_new.txt',
                "HGF Integrated PWPE2": base_folder+'HGF1730 PWPE2_mod_baked_new.txt',
                "HGF Integrated PWPE3": base_folder+'HGF1730 PWPE3_mod_baked_new.txt'
                }

for file in file_dict:
    pred = get_data(file_dict[file])
    try:
        pred = process_file(pred)   
    except:
        pred = process_file2(pred)
    #print(pred)
    plt.plot(pred,linestyle="",marker="o")
    plt.title(file)
    plt.legend(file)
    plt.show()

