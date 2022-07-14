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
import matplotlib.font_manager as font_manager
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


from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy import stats




# Attributes to (optionally) delete from each individual participant when loading the whole experiment
attrs_to_delete = [
                    'epochs',
                    'epochs_sound',
                    'epochs_ransac',
                   #'epochs_ransac_autoreject'
                   #'epochs_ransac_autoreject_unstandardised',
                   'events_copy','times','end_time', 
                   
                   #'surprise', 'surprise_category','surprise_predictions',
                   'sequence', 'tone_sequence',
                   'ecg_epochs', 'eog_epochs', 
                    #'ecg_events', 'ecg_indices', 'ecg_scores', 'eog_events', 'eog_indices', 'eog_scores', 
                    #'ransac',
                    'ar',

                   #'ica',
                   'ica_copy',
                   'raw','raw_filt',
                   'raw_cleaned',
                   #'difference_df'
                   ]


class Experiment():
    '''Full class for comparing participants' results'''
    
    def __init__(self, participant_strings=[], experiment_numbers=1, participant_pickles = None, condition=None):
        '''Initialise
        '''
        self.participants = participant_strings
        self.participants_data = participant_pickles # Optionally pump in a list of .pickles for each participant
        self.experiment_numbers = experiment_numbers
        self.condition_to_compare = condition
        
        # Needed to get metrics about the participants:
        self.all_ecg_excludes_ct= []
        self.all_eog_excludes_ct= []
        self.all_ecg_excludes   = []
        self.all_eog_excludes   = []
        self.all_reject_logs    = []
        self.all_reject_indices = []
        self.all_noisy          = []
        self.all_bad            = []
        self.all_flat           = []
        
        self.font = font_manager.FontProperties(family='Times New Roman',
                                               style='normal', size=FONTSIZE_LEGEND)


    def add_participants(self,participant_strings, is_adult_folder):
        '''Add in participants from their ID strings
        '''
        for participant_string in participant_strings:
            ptcp = Participant(is_adult_folder, p_id=participant_string, init=None, experiment_number=1)
            ptcp.basic_cleaning()
            ptcp.initiate_events()
            ptcp.more_processing()
            self.participants.append(ptcp)
        
    def in_memory(self, p_id):
        '''
        Input: a participant ID
        Output: whether that ID is in memory
        '''
        found_in_mem = False
        for participant in self.participants:
            if p_id in participant.p_id:
                found_in_mem = True         
        return found_in_mem        
        
    def add_participants_from_memory(self,ptcp_objects):
        '''
        Input: a list of Participant objects
        Adds these to memory if they aren't already in
        '''
        for participant in ptcp_objects:
            if not self.in_memory(participant.p_id):
                self.participants.append(participant)        
        
    def add_participants_from_disk(self,pickle_files):
        '''
        Input: pickle file paths
        Adds these to the experiment
        '''
        for pickle_file in pickle_files:
            with open(pickle_file, "rb") as f:
                participant = pickle.load(f)
            if not self.in_memory(participant.p_id):
                self.participants.append(participant)

        del participant
        gc.collect()


    def update_tracking(self, ptcp):
        '''Update tracking of processing results
        This is run after each participant is processed to minimise memory usage
        '''
        
        try: # Automated bad channels                           
            self.all_noisy.append(ptcp.auto_noisy_chs)
            self.all_flat.append(ptcp.auto_flat_chs)
        except:
            print("Missing bads info 1")
        try: # Bad channels
            self.all_bad.append(ptcp.epochs_ransac_autoreject.info['bads'])
        except:
            print("Missing bads info 2")            
        try: # Autoreject metrics
            self.all_reject_logs.append(ptcp.reject_log)
            self.all_reject_indices.append(ptcp.reject_indices)
            #self.epochs_ransac_autoreject.plot_drop_log()
        except:
            print("Missing AR or Ransac info")
        try: # ICA Metrics                 
            self.all_ecg_excludes_ct.append(ptcp.ecg_bad_ct)
            self.all_eog_excludes_ct.append(ptcp.eog_bad_ct)
            self.all_ecg_excludes.append(ptcp.ecg_bads)
            self.all_eog_excludes.append(ptcp.eog_bads)
        except:
            print("Missing ICA info")

        
    def match_channels(self):
        '''
        Load file that maps the adult sensors to appropriate child sensors 
        So we leave the child channels as 1,2,3,...125 and then we pick the adult channels that best match these
        This is necessary if e.g. plotting the topography of the brain signal as MEG001 in the adult system does not correspond to MEG001 in the child system (for example) in terms of where it is located in the brain

        It's not enough to just drop the non-comparable channels; we also have to fool MNE into thinking that the channels are comparable when e.g. doing a grand average as otherwise it'll drop the channels
        To do this we need to give the corresponding adult channels the same names as the children's (i.e. drop 126 to 160)
        < does MNE plot_topog know the location is different between systems ? >
        '''

        adult_sensors = scipy.io.loadmat(r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Lance\MEG\Matlab MEG\\labels_systemmatch_adult_sensors.mat')
        child_sensors = scipy.io.loadmat(r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Lance\MEG\Matlab MEG\\labels_systemmatch_child_sensors.mat')

        adult_sensors = list(adult_sensors.values())[3]
        child_sensors = list(child_sensors.values())[3]
        self.sensor_mapping_dct_simple = {}
        for r in range(0,len(adult_sensors)):
            item_adult = adult_sensors[r][0]
            item_child = child_sensors[r][0]
            num_adult = int(item_adult[0].replace("AG",""))
            num_child = int(item_child[0].replace("AG","")) 
            self.sensor_mapping_dct_simple["MEG "+str(num_adult).zfill(3)] = "MEG " + str(num_child).zfill(3)

    def align_channels(self):
        '''Perform a remapping of the adult channel names to child channel names to allow comparison of data across adult vs child machines
        '''
        self.match_channels()
        for ptcp in self.participants:
            #ptcp.align_channels(mapping=self.sensor_mapping_dct_simple)
            ptcp.align_channels_paul_individual()

    def find_age_groupings(self, system=None, num_age_groupings=int(100/PERCENTILES_AGE[1])):
        '''Creates num_age_groupings of participants based on age, each group equal sized
        If you specify the system, it only looks at participants scanned using that system ('Child' or 'Adult')
        '''
        if system==None:
            all_ages = [ptcp.age for ptcp in self.participants]
        else:
            all_ages = []
            for ptcp in self.participants:
                if (system=='Child' and ptcp.is_adult_system == False) or (system=='Adult' and ptcp.is_adult_system == True):          
                    all_ages.append(ptcp.age)
        ages_linear = np.arange(min(all_ages),max(all_ages),num_age_groupings)
        age_bounds = [(ages_linear[r], ages_linear[r+1]) for r in range(0,len(ages_linear)-1)] #  e.g.  [(2,5),(5,7),(6,8),(8,21),(21,60)]
        age_boundaries = equalObs(all_ages, num_age_groupings)
        self.age_bounds = [(float("{:.2f}".format(age_boundaries[a])),float("{:.2f}".format(age_boundaries[a+1]))) for a in range(0,num_age_groupings)] 
        print("Age groupings are ", self.age_bounds)




    def analyse_age_groupings(self,gfp_mode=True, first_row=True):
        '''Analyse condition A vs condition B in each age group
        If gfp mode, analyses gfps, else analyses ERFs
        
        # When plotting bottom row, set first_row=False:
        i) Adds back time (ms) as x axis label
        ii) Removes legend and title
        Changes params min_plot, max_plot with this
        '''

        import matplotlib
        from matplotlib.ticker import MaxNLocator

        first_row = False # True
        min_plot = 130 # 0 # crop_start*1000
        max_plot = 190 # 300 #crop_end*1000
        yminplot_1 = 0
        ymaxplot_1 = 25
        y_ticks = [0, 5, 10, 15, 20, 25]
        y_ticks2 = [-8,-4, 0, 4, 8]
        yminplot_2 = -7
        ymaxplot_2 = 7


        if len(self.participants) > 6:
            self.find_age_groupings()

            # Condition in various age group ranges
            evoked_all_A = {}
            evoked_all_B = {}
            evoked_everything = {}            
            for age_min, age_max in self.age_bounds:
                evoked_all_A[(age_min,age_max)] = []
                evoked_all_B[(age_min,age_max)] = []
                evoked_everything[(age_min,age_max)] = []                
                for ptcp in self.participants:
                    if ptcp.age >= age_min and ptcp.age < age_max:
                        #if len(ptcp.evoked_generic[ptcp.cond_A].to_data_frame()) == NUM_EPOCH_SAMPLES and len(ptcp.evoked_generic[ptcp.cond_B].to_data_frame()) == NUM_EPOCH_SAMPLES:
                        evoked_all_A[(age_min,age_max)].append(ptcp.evoked_generic[ptcp.cond_A])
                        evoked_all_B[(age_min,age_max)].append(ptcp.evoked_generic[ptcp.cond_B])
                        evoked_everything[(age_min,age_max)].append(ptcp.evoked_generic[ptcp.cond_A])
                        evoked_everything[(age_min,age_max)].append(ptcp.evoked_generic[ptcp.cond_B])                              

            for key in evoked_all_A.keys():
                print("Age grouping ", key,len(evoked_all_A[key]))
            self.grand_averages_A = {}        
            self.grand_averages_B = {}        
            for key in evoked_all_A:
                self.grand_averages_A[key] = mne.grand_average(evoked_all_A[key])
                self.grand_averages_A_df = self.grand_averages_A[key].to_data_frame()
                self.grand_averages_A_gfp_vs_time = RMS_DF(self.grand_averages_A_df)   
                latency = find_latency(self.grand_averages_A_df, self.grand_averages_A_gfp_vs_time )
                print(key,latency)                                
            for key in evoked_all_B:
                self.grand_averages_B[key]  = mne.grand_average(evoked_all_B[key])       
                self.grand_averages_B_df = self.grand_averages_B[key].to_data_frame()
                self.grand_averages_B_gfp_vs_time = RMS_DF(self.grand_averages_B_df)    
                latency = find_latency(self.grand_averages_B_df, self.grand_averages_B_gfp_vs_time )
                print(key,latency)

            self.latencies_by_age = {}
            for age_min, age_max in self.age_bounds:
                latencies = []
                self.latencies_by_age[(age_min, age_max)] = []
                for ptcp in evoked_everything[(age_min,age_max)]:
                    ptcp_df = ptcp.to_data_frame()
                    ptcp_df_gfp = RMS_DF(ptcp_df)   
                    latency = find_latency(ptcp_df, ptcp_df_gfp)
                    self.latencies_by_age[(age_min, age_max)].append(latency)
                print((age_min, age_max), np.mean(self.latencies_by_age[(age_min, age_max)]),np.std(self.latencies_by_age[(age_min, age_max)]))


            from pingouin import ttest


            print(ttest(latencies_by_age[(3.2, 5.1)], latencies_by_age[(5.1, 6.8)], correction = True))
            for r in range(0,len(latencies_by_age.keys())-1):
                key_1, key_2 = latencies_by_age.keys()[r], atencies_by_age.keys()[r+1]
                stats_tests = ttest(latencies_by_age[key_1], latencies_by_age[(key_2)], correction = True).to_dict()
                stats_df = {}
                for stat in stats_tests.keys():
                    try:
                        if type(float(stats_tests[stat]['T-test'])) == type(float(0.0034)):
                            stats_df[stat] = stats_tests[stat]['T-test']
                            #print(stats_df)
                    except:
                        pass
                stats_df = pd.DataFrame(list(stats_df.items()))
                stats_df



            #NEW
            fig=plt.figure()
            ax=fig.gca() # add_subplot(111)

            # Compare condition A for each age grouping       
            plt_legend = []
            for key in self.grand_averages_A.keys():
                #print("Age grouping ", key)
                if gfp_mode:
  
                    
                    ax.plot(self.grand_averages_A_df['time'].values, self.grand_averages_A_gfp_vs_time)
                    #plt.plot(self.grand_averages_A_df['time'].values, self.grand_averages_A_gfp_vs_time)
                    #ax.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS, xmax=EPOCH_END_ANALYSIS*SEC_TO_MS)
                    ax.set_xlim([PLOT_START*SEC_TO_MS,PLOT_END*SEC_TO_MS])
                    #ax.set_xlim(xmin=min_plot,xmax=max_plot)
                    ax.set_ylim(ymin=yminplot_1, ymax=ymaxplot_1)    # HARDCODE
                    
                    ax.set_yticks(y_ticks)
                    ax.set_yticklabels(y_ticks)
                    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
                    ax.yaxis.set_major_locator(MaxNLocator(len(y_ticks))) 

                    # Hide the right and top spines
                    #ax.spines['right'].set_visible(False)
                    #ax.spines['top'].set_visible(False)
                    # Only show ticks on the left and bottom spines
                    ax.yaxis.set_ticks_position('left')
                    ax.xaxis.set_ticks_position('bottom')
                    #ax.set_yticklabels([])
                    matplotlib.rcParams["axes.spines.right"] = False
                    matplotlib.rcParams["axes.spines.top"] = False


                else: # Meaningless plot
                    #Somewhat meaningless 'sum' column
                    self.grand_averages_A_df = sum_df(self.grand_averages_A_df,mode='plain')
                    # Plot ERFs
                    #plt.plot(self.grand_averages_A_df['time'].values, self.grand_averages_A_df['sum'].values)
                    #plt.legend(["ERF Condition:" +str(self.participants[0].cond_A), "ERF Condition:"+str(self.participants[0].cond_B)], bbox_to_anchor=(0.75, 1.15), ncol=2, fontsize=FONTSIZE_LEGEND)
                #if first_row:
                #plt_legend.append("Age Range "+str(key))

            if gfp_mode:
                if first_row == True:
                    pass
                    #plt.title("GFP for low surprise trials", fontsize=FONTSIZE_TITLE)
                plt.ylabel("GFP (fT)", fontsize=FONTSIZE_LABELS)      
                #if not first_row:
                plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)                
                #if first_row:
                plt.legend(plt_legend, fontsize=FONTSIZE_LEGEND, loc='best', bbox_to_anchor=(0.6, 1.2))
                plt.savefig("GFPs Age Groups Cond A "+str(int(EPOCH_START_ANALYSIS*SEC_TO_MS))+"-"+str(int(EPOCH_END_ANALYSIS*SEC_TO_MS))+" ms.jpg" , bbox_inches='tight', dpi=DPI) # 
                plt.show()                              
            else:
                pass
                # plt.title("ERFs Condition A", fontsize=FONTSIZE_TITLE)
                # plt.ylabel("ERF (fT)", fontsize=FONTSIZE_LABELS)  
                #plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)                
                # plt.legend(plt_legend, fontsize=FONTSIZE_LEGEND, loc='best', bbox_to_anchor=(0.6, 1.2)))
                # plt.show()     
            

            #NEW
            fig=plt.figure()
            ax=fig.gca() # ax=fig.add_subplot(111)
            # Compare condition B for each age grouping
            plt_legend = []
            for key in self.grand_averages_B.keys():
                print("Age grouping ", key)                
                self.grand_averages_B_df = self.grand_averages_B[key].to_data_frame()

                if gfp_mode:
                    self.grand_averages_B_gfp_vs_time = RMS_DF(self.grand_averages_B_df)    
                    #plt.plot(self.grand_averages_B_df['time'].values,  self.grand_averages_B_gfp_vs_time)
                    #plt.legend(["GFP Condition:" +str(self.participants[0].cond_A), "GFP Condition:"+str(self.participants[0].cond_B)], bbox_to_anchor=(0.75, 1.15), ncol=2, fontsize=FONTSIZE_LEGEND)                    
                    
                    #plt.plot(self.grand_averages_B_df['time'].values,self.grand_averages_B_gfp_vs_time)
                    ax.plot(self.grand_averages_B_df['time'].values, self.grand_averages_B_gfp_vs_time)
                    #ax.set_xlim(xmin=min_plot,xmax=max_plot)
                    ax.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS,xmax=EPOCH_END_ANALYSIS*SEC_TO_MS) 
                    ax.set_xlim(xmin=PLOT_START*SEC_TO_MS,xmax=PLOT_END*SEC_TO_MS) 
                    ax.set_ylim(ymin=yminplot_1, ymax=ymaxplot_1)    # HARDCODE 
                    ax.set_yticks(y_ticks) 
                    ax.set_yticklabels(y_ticks)                
                    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
                    ax.yaxis.set_major_locator(MaxNLocator(len(y_ticks))) 

                    # Hide the right and top spines
                    #ax.spines['right'].set_visible(False)
                    #ax.spines['top'].set_visible(False)
                    # Only show ticks on the left and bottom spines
                    ax.yaxis.set_ticks_position('left')
                    ax.xaxis.set_ticks_position('bottom')
                    ax.set_yticklabels([])
                    matplotlib.rcParams["axes.spines.right"] = False
                    matplotlib.rcParams["axes.spines.top"] = False

                else: # Meaningless plot
                    #Somewhat meaningless 'sum' column
                    self.grand_averages_B_df = sum_df(self.grand_averages_B_df,mode='plain')
                    # Plot ERFs
                    #plt.plot(self.grand_averages_B_df['time'].values, self.grand_averages_B_df['sum'].values)
                    #plt.legend(["ERF Condition:" +str(self.participants[0].cond_A), "ERF Condition:"+str(self.participants[0].cond_B)], bbox_to_anchor=(0.75, 1.15), ncol=2, fontsize=FONTSIZE_LEGEND)
                #plt_legend.append("Age Range "+str(key))
                
                
            if gfp_mode:
                if first_row == True:
                    pass
                    #plt.title("GFP for high surprise trials", fontsize=FONTSIZE_TITLE)
                #plt.ylabel("GFP (fT)", fontsize=FONTSIZE_LABELS)   # Hide label as it's obvious from the other graph 
                #if not first_row:
                plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)                
                #plt.legend(plt_legend, fontsize=FONTSIZE_LEGEND)
                plt.savefig("GFPs Age Groups Cond B "+str(int(EPOCH_START_ANALYSIS*SEC_TO_MS))+"-"+str(int(EPOCH_END_ANALYSIS*SEC_TO_MS))+" ms.jpg" , bbox_inches='tight', dpi=DPI) # 
                plt.show()                                  
            else:
                pass
                # plt.title("ERFs for high surprise trials", fontsize=FONTSIZE_TITLE)
                # plt.ylabel("ERF (fT)", fontsize=FONTSIZE_LABELS)  
                # plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)                
                # plt.legend(plt_legend, fontsize=FONTSIZE_LEGEND)
                # plt.show()    
            

            #NEW
            fig=plt.figure()
            ax=fig.gca()  # ax=fig.add_subplot(111)


              

            # Compare difference in condition B to condition A for each age grouping
            plt_legend = []
            for key in self.grand_averages_A.keys():
                print("Age grouping ", key)
                if gfp_mode:   # Option A - # Plot difference in GFPs
                    self.grand_averages_A_df = self.grand_averages_A[key].to_data_frame()
                    self.grand_averages_B_df = self.grand_averages_B[key].to_data_frame()
                    self.grand_averages_A_gfp_vs_time = RMS_DF(self.grand_averages_A_df)
                    self.grand_averages_B_gfp_vs_time = RMS_DF(self.grand_averages_B_df)
                    gfp_diff = [self.grand_averages_B_gfp_vs_time[r] - self.grand_averages_A_gfp_vs_time[r] for r in range(0,len(self.grand_averages_B_gfp_vs_time))]

                    #plt.plot(self.grand_averages_A_df['time'].values,gfp_diff)
                    ax.plot(self.grand_averages_A_df['time'].values, gfp_diff)
                    ax.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS,xmax=EPOCH_END_ANALYSIS*SEC_TO_MS) 
                    ax.set_xlim(xmin=PLOT_START*SEC_TO_MS,xmax=PLOT_END*SEC_TO_MS) 
                    #ax.set_xlim(xmin=min_plot,xmax=max_plot)
                    ax.set_ylim(ymin=yminplot_2, ymax=ymaxplot_2)    # HARDCODE
                    ax.set_yticks(y_ticks2)
                    ax.set_yticklabels(y_ticks2)                
                    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
                    ax.yaxis.set_major_locator(MaxNLocator(len(y_ticks2))) 

                    # Hide the right and top spines
                    #ax.spines['right'].set_visible(False)
                    #ax.spines['top'].set_visible(False)
                    # Only show ticks on the left and bottom spines
                    ax.yaxis.set_ticks_position('left')
                    ax.xaxis.set_ticks_position('bottom')
                    matplotlib.rcParams["axes.spines.right"] = False
                    matplotlib.rcParams["axes.spines.top"] = False

                else:           # Option B - Plot difference in ERFs (MMR)
                    #Somewhat meaningless 'sum' column
                    self.grand_averages_A_df = self.grand_averages_A[key].to_data_frame()
                    self.grand_averages_B_df = self.grand_averages_B[key].to_data_frame()                    
                    self.grand_averages_A_df = sum_df(self.grand_averages_A_df,mode='plain')
                    self.grand_averages_B_df = sum_df(self.grand_averages_B_df,mode='plain')

                self.grand_average_cond_AB_diff = self.grand_averages_B_df - self.grand_averages_A_df
                self.grand_average_cond_AB_diff['time'] = self.grand_averages_A_df['time']
                #plt.plot(self.grand_average_cond_AB_diff['time'].values,self.grand_average_cond_AB_diff['sum'].values) 
                #plt_legend.append("Age Range "+str(key))

            if gfp_mode:
                if first_row == True:
                    pass
                    #plt.title("GFP for high surprise trials\n"+"minus GFP for low surprise trials", fontsize=FONTSIZE_TITLE)
                plt.ylabel("GFP Difference (fT)", fontsize=FONTSIZE_LABELS)                  
            else:
                pass
                # plt.title("MMR (sum) Condition B minus condition A", fontsize=FONTSIZE_TITLE)
                # plt.ylabel("ERF Difference (fT)", fontsize=FONTSIZE_LABELS)  
            
            #if not first_row:
            plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
            #plt.legend(plt_legend, fontsize=FONTSIZE_LEGEND, loc='best', bbox_to_anchor=(1, 0.5))
            plt.savefig("GFPs Age Groups (B-A) "+str(int(EPOCH_START_ANALYSIS*SEC_TO_MS))+"-"+str(int(EPOCH_END_ANALYSIS*SEC_TO_MS))+" ms.jpg" , bbox_inches='tight', dpi=DPI) # 
            plt.show() 

            del evoked_all_A, evoked_all_B
            gc.collect()
        else:
            print("You only have one participant, cannot do groupings")


    def group_participants_on_age(self,num_per_group=None, absolute=False, system=None, age_bounds_low=None, age_bounds_high = None, age_cutoff=AGE_CUTOFF):
        '''Splits participants into two groups: young and old
        Makes self.group_A (youngest N persons), self.group_B (oldest N persons)
        Uses half split if num_per_group == None
        
        if absolute == True, divides them not into equally-sized groups but into potentially unequally-sized groups based on absolute age criterion
        Not fully tested as to it not causing errors later on

        if system='Child' or 'Adult', only picks participants scanned in the child/adult system respectively

        e.g. absolute=True, system='Child',age_cutoff = 6
        -> self.youngest = people < 6 scanned in child system
        -> self.oldest = people >= 6 scanned in child system
        '''
        if absolute == True:
            self.youngest = []
            self.oldest = []
            for participant in self.participants:
                if age_bounds_low==None and age_bounds_high == None: # Then the young and old age groups respectively are just younger or older than the cutoff
                    if participant.age <= age_cutoff:
                        if system==None or (system=='Child' and ptcp.is_adult_system == False) or (system=='Adult' and ptcp.is_adult_system == True):    
                            self.youngest.append(participant.p_id)
                    else:
                        if system==None or (system=='Child' and ptcp.is_adult_system == False) or (system=='Adult' and ptcp.is_adult_system == True):                    
                            self.oldest.append(participant.p_id)
                else:


                    if age_bounds_low[1] > age_bounds_high[0]:
                        self.warn("Age bounds don't make sense, should be mutually exclusive")
                    if participant.age >= age_bounds_low[0] and participant.age <= age_bounds_low[1]:
                        if system==None or (system=='Child' and ptcp.is_adult_system == False) or (system=='Adult' and ptcp.is_adult_system == True):    
                            self.youngest.append(participant.p_id)
                    elif participant.age > age_bounds_high[0] and participant.age <= age_bounds_high[1]:
                        if system==None or (system=='Child' and ptcp.is_adult_system == False) or (system=='Adult' and ptcp.is_adult_system == True):                    
                            self.oldest.append(participant.p_id)

            print("Young: ", self.youngest)
            print("Old: ", self.oldest)
            self.group_participants(self.youngest,self.oldest)

        else: # Then define young and old relatively
            if not num_per_group:
                num_per_group = int(len(self.participants)*0.5)
            num_per_group = min(int(num_per_group),int(len(self.participants)*0.5)) # Logical value check

            # Method A : memory-efficient
            self.participant_age_dict = {}
            for participant in self.participants:
                if participant.age!= None:
                    self.participant_age_dict[participant.p_id] = participant.age
                else:
                    print("%%%%%% Unknown age for participant %%%%%%: ", participant.p_id )
            self.participant_strings_by_age = {k: v for k, v in sorted(self.participant_age_dict.items(), key=lambda item: item[1])}

            # To sort the list in place...
            self.participants.sort(key=lambda x: x.age)#, reverse=True)
            youngest_half, oldest_half = split_list(list(self.participant_strings_by_age.keys())) 
            print("Youngest half ", youngest_half, "Oldest half ", oldest_half)
            self.youngest = youngest_half[0:num_per_group]
            self.oldest = oldest_half[-num_per_group:]
            print("Youngest: ", self.youngest)
            print("Oldest: ", self.oldest)
            self.group_participants(self.youngest,self.oldest)
            self.record_age_bounds()

    def group_participants(self,group_A_ids,group_B_ids):
        '''Splits participants into two groups'''
        self.group_A = group_A_ids
        self.group_B = group_B_ids

    def record_age_bounds(self):
        self.age_young_min = 999999
        self.age_young_max = -999999
        self.age_old_min = 999999
        self.age_old_max = -999999
        for participant in self.participants:
            if participant.p_id in self.youngest:
                if participant.age < self.age_young_min: 
                    self.age_young_min = participant.age
                if participant.age > self.age_young_max: 
                    self.age_young_max = participant.age
            if participant.p_id in self.oldest:
                if participant.age < self.age_old_min: 
                    self.age_old_min = participant.age
                if participant.age > self.age_old_max: 
                    self.age_old_max = participant.age                    



        
    def compare_ERFs_conditions(self,cond_A,cond_B):
        ''' Computes ERFs and compares for two specific conditions over ALL participants (i.e. grand average ERFs)'''
        self.cond_A_evoked = []
        self.cond_B_evoked = []
        r = 0
        for ptcp in self.participants:
            if r % 10 == 0:
                print("Compare ERFs loop: participant #: ", r)
            c_a = ptcp.evoked_generic[cond_A]
            c_b = ptcp.evoked_generic[cond_B]
            #Somewhat meaningless 'sum' column
            c_a_df = sum_df(c_a.to_data_frame())
            c_b_df = sum_df(c_b.to_data_frame())
            sum_a = c_a_df['sum'].values
            sum_b = c_b_df['sum'].values
            if sum([1 for x in sum_a if math.isnan(x)]) >0 or sum([1 for x in sum_b if math.isnan(x)]) >0:
                self.warn("NANs present for " + str(self.p_id) +" not adding")
            else:
                self.cond_A_evoked.append(c_a) 
                self.cond_B_evoked.append(c_b)
            for var_to_delete in ['c_a','c_b']: #",'x','y']:
                globals().pop(var_to_delete, None)
            gc.collect()
            r+=1

        self.grand_average_cond_A = mne.grand_average(self.cond_A_evoked)
        clear_output()
        self.grand_average_cond_B = mne.grand_average(self.cond_B_evoked)
        clear_output()
        self.grand_average_cond_A_df = self.grand_average_cond_A.to_data_frame()
        self.grand_average_cond_B_df = self.grand_average_cond_B.to_data_frame()
        #Somewhat meaningless 'sum' column
        self.grand_average_cond_A_df = sum_df(self.grand_average_cond_A_df) #, bads=self.bads)
        self.grand_average_cond_B_df = sum_df(self.grand_average_cond_B_df) #, bads=self.bads)   

                  
        clear_output()
        # Plot ERFs
        # plt.plot(self.grand_average_cond_A_df['time'].values, self.grand_average_cond_A_df['sum'].values)
        # plt.plot(self.grand_average_cond_B_df['time'].values, self.grand_average_cond_B_df['sum'].values)
        # plt.legend(["ERF low surprise condition:" +str(cond_A), "ERF high surprise condition:"+str(cond_B)], bbox_to_anchor=(0.75, 1.15), ncol=2, fontsize=FONTSIZE_LEGEND)
        # plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
        # plt.ylabel("Sum of ERFs across channels (fT)", fontsize=FONTSIZE_LABELS)
        # plt.title("Sum of ERFs across channels, in low (blue) & high (orange) expected surprise conditions", fontsize=FONTSIZE_TITLE)
        # plt.show()

        # Plot difference in ERFs (MMR)
        self.grand_average_cond_AB_diff = self.grand_average_cond_B_df - self.grand_average_cond_A_df
        self.grand_average_cond_AB_diff['time'] = self.grand_average_cond_A_df['time']
        # plt.plot(self.grand_average_cond_AB_diff['time'].values, self.grand_average_cond_AB_diff['sum'].values)
        # plt.legend(["ERF sum, high surprise condition:" +str(cond_B) +" *minus* low surprise condition:"+str(cond_A)], bbox_to_anchor=(0.75, 1.15), ncol=2, fontsize=FONTSIZE_LEGEND)
        # plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
        # plt.ylabel("Difference of sum of ERFs across channels (fT)", fontsize=FONTSIZE_LABELS)
        # plt.title("MMR Cond B (high E[surprise]) - A (low E[surprise])", fontsize=FONTSIZE_TITLE)
        # plt.show()    
        
        
        
    def average_evoked(self,participant_strings,condition=None):
        '''Averages the evoked for condition B (or A)   (for participants with p_ids in participant_strings)
        '''

        evokeds = []
        r = 0
        for ptcp in self.participants:
            if r % 15 == 0: # Progress update
                print("Average evoked loop: participant #: ", r)
            if str(ptcp.p_id) in participant_strings:
                conds_to_compare = ptcp.conds_to_compare
                events_to_tag = ptcp.events_to_tag

                if condition==None:
                    #if len(ptcp.evoked_all.to_data_frame()) == NUM_EPOCH_SAMPLES: # Enforce that the evokeds are the same length else we cannot compare conditions    
                    evoked_object = ptcp.evoked_all
                    evokeds.append(evoked_object) # Compare ERFs for a specific condition (condition B) in Group A to Group B
                elif condition=='B':
                    #if len(ptcp.evoked_generic[self.condition_to_compare].to_data_frame()) == NUM_EPOCH_SAMPLES: # Enforce that the evokeds are the same length else we cannot compare conditions    
                    evoked_object = ptcp.evoked_generic[conds_to_compare[events_to_tag[1]][1]] # self.condition_to_compare]  
                    evokeds.append(evoked_object) # Compare ERFs for a specific condition (condition B) in Group A to Group B
                elif condition=='A':
                    #if len(ptcp.evoked_generic[conds_to_compare[events_to_tag[1]][0]].to_data_frame()) == NUM_EPOCH_SAMPLES: # Enforce that the evokeds are the same length else we cannot compare conditions  
                    evoked_object = ptcp.evoked_generic[conds_to_compare[events_to_tag[1]][0]] 
                    evokeds.append(evoked_object) # Compare ERFs for a specific condition  (condition A) in Group A to Group B 
            r+=1

        if len(evokeds) > 0:
            grand_average = mne.grand_average(evokeds)
            del evoked_object, evokeds
            gc.collect()
            return grand_average
        else:
            print("Problem getting evokeds")


    def compare_group(self, condition=None):
        '''Compares Group 1 ERF (for condition B) minus Group 2 ERF (for condition B)'''

        if not hasattr(self,'group_A_avg'):
            self.group_A_avg = {}
        if not hasattr(self,'group_B_avg'):
            self.group_B_avg = {}
        
        # Ensure cropping is done
        crop_start = EPOCH_START_ANALYSIS
        crop_end = EPOCH_END_ANALYSIS
        for ptcp in self.participants:
            for ptcp in self.participants:
                ptcp.epochs_ransac_autoreject.crop(tmin=crop_start, tmax=crop_end, include_tmax=True)
                ptcp.evoked_generic[ptcp.cond_A].crop(tmin=crop_start, tmax=crop_end, include_tmax=True)
                ptcp.evoked_generic[ptcp.cond_B].crop(tmin=crop_start, tmax=crop_end, include_tmax=True)
                ptcp.evoked_all.crop(tmin=crop_start, tmax=crop_end, include_tmax=True)   


        # Group 1 grand average
        self.group_A_avg[condition] = self.average_evoked(self.group_A,condition=condition) # Get the average for condition B amongst every member of Group A
        self.group_B_avg[condition] = self.average_evoked(self.group_B,condition=condition) # Get the average for condition B amongst every member of Group B
        self.group_A_avg_df = self.group_A_avg[condition].to_data_frame()
        self.group_B_avg_df = self.group_B_avg[condition].to_data_frame()

        #Somewhat meaningless 'sum' column
        self.group_A_avg_df = sum_df(self.group_A_avg_df,mode='plain')  
        self.group_B_avg_df = sum_df(self.group_B_avg_df,mode='plain')              

        #clear_output()
        # Plot ERFs
        # plt.plot(self.group_A_avg_df['time'].values, self.group_A_avg_df['sum'].values)
        # plt.plot(self.group_B_avg_df['time'].values, self.group_B_avg_df['sum'].values)
        # plt.legend(["ERF Group A (younger)", "ERF Group B (older) for condition"+str(condition)], bbox_to_anchor=(0.75, 1.15), ncol=2, fontsize=FONTSIZE_LEGEND)
        # plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
        # plt.ylabel("ERF (fT)", fontsize=FONTSIZE_LABELS) 
        # plt.title("ERF Group A (younger), B (older) for condition "+str(condition), fontsize=FONTSIZE_TITLE)
        # plt.show()     

        ## Calculate difference between groups
        # Method A
        self.group_diff_df = self.group_A_avg_df - self.group_B_avg_df
        self.group_diff_df['time'] = self.group_A_avg_df['time']

    def compare_group_mmrs(self):
        '''Compares Group 1 MMR (condition B ERF minus condition A ERF) minus Group 2 MMR (condition B ERF minus condition A ERF)'''

        # Group 1 grand average
        self.group_A_avg_condB = self.average_evoked(self.group_A,condition='B') # Get the average for condition B amongst every member of Group A        
        self.group_A_avg_condA = self.average_evoked(self.group_A,condition='A') # Get the average for condition A amongst every member of Group A
        self.group_B_avg_condB = self.average_evoked(self.group_B,condition='B') # Get the average for condition B amongst every member of Group B
        self.group_B_avg_condA = self.average_evoked(self.group_B,condition='A') # Get the average for condition A amongst every member of Group B

        # Convert to DFs
        self.group_A_avg_condB_df = self.group_A_avg_condB.to_data_frame()
        self.group_A_avg_condA_df = self.group_A_avg_condA.to_data_frame()
        self.group_B_avg_condB_df = self.group_B_avg_condB.to_data_frame()
        self.group_B_avg_condA_df = self.group_B_avg_condA.to_data_frame()


        self.group_A_avg_MMR_df = self.group_A_avg_condB_df - self.group_A_avg_condA_df
        self.group_B_avg_MMR_df = self.group_B_avg_condB_df - self.group_B_avg_condA_df        
        # Somewhat meaningless 'sum' column
        self.group_A_avg_MMR_df = sum_df(self.group_A_avg_MMR_df,mode='plain')  
        self.group_B_avg_MMR_df = sum_df(self.group_B_avg_MMR_df,mode='plain')              

        #clear_output()
        # Plot MMRs
        self.group_A_avg_MMR_df['time'] = self.group_A_avg_condB_df['time'].values
        self.group_B_avg_MMR_df['time'] = self.group_B_avg_condB_df['time'].values



        plt = plot_module(  xs = [self.group_A_avg_MMR_df['time'].values, self.group_B_avg_MMR_df['time'].values],
                            ys = [self.group_A_avg_MMR_df['sum'].values,  self.group_B_avg_MMR_df['sum'].values],
                            xlabel = "Time (ms)", ylabel="MMR (fT)",
                            legends=["MMR Group A (younger)", "MMR Group B (older)"],
                            title="MMRs for Group A (younger) and Group B (younger)",
                            font_object = self.font,
                            x_min_to_show = PLOT_START*SEC_TO_MS,
                            x_max_to_show = PLOT_END*SEC_TO_MS,                              
                            #x_min_to_show = EPOCH_START_ANALYSIS*SEC_TO_MS,
                            #x_max_to_show = EPOCH_END_ANALYSIS*SEC_TO_MS                            
                            )

        # plt.plot(self.group_A_avg_MMR_df['time'].values, self.group_A_avg_MMR_df['sum'].values)
        # plt.plot(self.group_B_avg_MMR_df['time'].values, self.group_B_avg_MMR_df['sum'].values)
        # #plt.xlim([min(self.group_A_avg_MMR_df['time'].values),max(self.group_A_avg_MMR_df['time'].values)])
        # plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed
        # plt.legend(["MMR Group A (younger)", "MMR Group B (older)"], bbox_to_anchor=(0.75, 1.15), ncol=2, fontsize=FONTSIZE_LEGEND)
        # plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
        # plt.ylabel("MMR (fT)", fontsize=FONTSIZE_LABELS) 
        # plt.title("MMRs for Group A (younger) and Group B (younger)", fontsize=FONTSIZE_TITLE)
        plt.show()     

        ## Calculate difference in the MMR between groups
        # Method A
        self.group_diff_MMR_df = self.group_A_avg_MMR_df - self.group_B_avg_MMR_df
        self.group_diff_MMR_df['time'] = self.group_A_avg_MMR_df['time']




        # Plot difference in the MMR between groups

        plt = plot_module(  xs = [self.group_diff_MMR_df['time'].values],
                            ys = [self.group_diff_MMR_df['sum'].values],
                            xlabel = "Time (ms)", ylabel="MMR (fT)",
                            legends=["MMR difference (fT)"],
                            title="MMR Group A (younger) minus MMR Group B (older)",
                            font_object = self.font,
                            x_min_to_show = PLOT_START*SEC_TO_MS,
                            x_max_to_show = PLOT_END*SEC_TO_MS,                              
                            #x_min_to_show = EPOCH_START_ANALYSIS*SEC_TO_MS,
                            #x_max_to_show = EPOCH_END_ANALYSIS*SEC_TO_MS                            
                            )

        # plt.plot(self.group_diff_MMR_df['time'].values, self.group_diff_MMR_df['sum'].values)
        # #plt.xlim([min(self.group_diff_MMR_df['time'].values),max(self.group_diff_MMR_df['time'].values)])
        # plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed        
        # plt.legend(["Group A MMR\n"+str("minus Group B MMR")], bbox_to_anchor=(0.75, 1.15), ncol=2, fontsize=FONTSIZE_LEGEND)
        # plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
        # plt.ylabel("MMR difference (fT)", fontsize=FONTSIZE_LABELS)
        # plt.title("MMR Group A (younger) minus MMR Group B (older)", fontsize=FONTSIZE_TITLE)
        plt.show() 

     
    def remove_unwanted_data(self):
        
        '''
        Remove 
        a) select participants ('remove_ids')
        b) participants whose sound delay was not determined accurately (remove_bad_sound_delay)
        c) participants if they have a different condition to others in the group remove_wrong_condition/remove_missing_condition )
        c) channels from participants which we can't analyse (remove_missing_channels)
        or participants with errors  (remove_wrong_timespan)
        '''
        self.removed = {}
        
        self.remove_ids()
        self.remove_bad_sound_delay()
        self.rerun_if_wrong_condition()
        self.remove_wrong_condition() # Delete any which failed to re-run
        
        self.remove_missing_condition()
        self.remove_missing_channels()
        self.remove_surplus_channels()
        #self.remove_wrong_timespan()

    def remove_ids(self,ids=['9000','9001','9011']):
        '''Remove participants based on p_id
        '''    
        r = 0
        indices_to_delete = []
        for ptcp in self.participants:
            if ptcp.p_id in ids:
                indices_to_delete.append(r)
                print("Remove IDs loop: ptcp ", r, ptcp.p_id, ptcp.cond_A,ptcp.cond_B)
            r+=1
        indices_to_delete = list(set(indices_to_delete))
        print("Will delete based on p_id ", indices_to_delete)
        self.removed['IDs'] = indices_to_delete
        for index in sorted(indices_to_delete, reverse=True):
            del self.participants[index]

    def remove_bad_sound_delay(self):
        '''Remove participants whose median sound delay is nonsensical'''
        r = 0
        indices_to_delete = []
        for ptcp in self.participants:
            if ptcp.sound_delay > 0.3 or ptcp.sound_delay < 0:
                indices_to_delete.append(r)
                print("Remove based on sound delay loop: ptcp #:", r, "ID: ", ptcp.p_id, "Sound delay sec", ptcp.sound_delay)
            r+=1
        indices_to_delete = list(set(indices_to_delete))
        print("Will delete based on sound delay ", indices_to_delete)
        self.removed['SoundDelays'] = indices_to_delete
        for index in sorted(indices_to_delete, reverse=True):
            del self.participants[index]        

    def rerun_if_wrong_condition(self):
        '''Re-run participants whose pipeline was run analysing a different condition
        '''
        try:
            indices_to_rerun = []
            for ptcp in self.participants:
                if ptcp.cond_A != conds_to_compare[events_to_tag[1]][0] or ptcp.cond_B != conds_to_compare[events_to_tag[1]][1]: #or sorted(list(ptcp.evoked_generic.keys())) != [conds_to_compare[events_to_tag[1]][0],conds_to_compare[events_to_tag[1]][1]]:
                    ptcp.rerun_with_new_conditions(events_to_tag_dct = events_to_tag_rerun, conds_to_compare_dct=conds_to_compare)
                    print("Reran due to wrong condition ", ptcp.p_id, "New condition A: ", ptcp.cond_A, "New condition B: ", ptcp.cond_B, conds_to_compare[events_to_tag[1]][0], conds_to_compare[events_to_tag[1]][1])
        except Exception as e:
            print("Failed to re-run, error: ", e)
           
    def remove_wrong_condition(self):
        '''Remove participants whose pipeline was run analysing a different condition, and which got an error on attempt to re-run
        '''
        r = 0
        indices_to_delete = []
        for ptcp in self.participants:
            if ptcp.cond_A != conds_to_compare[events_to_tag[1]][0] or ptcp.cond_B != conds_to_compare[events_to_tag[1]][1]: # or sorted(list(ptcp.evoked_generic.keys())) != [conds_to_compare[events_to_tag[1]][0],conds_to_compare[events_to_tag[1]][1]]:
                indices_to_delete.append(r)
                print("Remove wrong condition loop: ptcp # ", r, " ID : ", ptcp.p_id, "Condition A: ", ptcp.cond_A, "Condition B: ", ptcp.cond_B, conds_to_compare[events_to_tag[1]][0], conds_to_compare[events_to_tag[1]][1])
            r+=1
        indices_to_delete = list(set(indices_to_delete))
        indices_to_delete = sorted(indices_to_delete, reverse=True)  
        self.removed['Wrong Condition'] = indices_to_delete
        print("Will delete based on wrong conditions", indices_to_delete)
        for index in indices_to_delete:
            del self.participants[index]
            
    def get_common_channels(self,evoked_A,evoked_B):
        ''' Identify which channels are in common between two evoked objects
        Used to allow comparison between those scanned in adult and child systems
        '''
        chs_A = evoked_A.info['ch_names']
        chs_B = evoked_B.info['ch_names']
        not_common = [x for x in chs_A if x not in chs_B]
        not_common = not_common + [x for x in chs_B if x not in chs_A]
        return not_common
    
    def identify_missing_channels(self):
        '''Identify any equivalents of child channels missing in any participant's data. If missing in one can't do the group analysis using it.
        '''
        ch_names = []
        for ch in range(1,num_meg_chs_child):
            ch_name = "MEG "+str(ch).zfill(3)
            ch_names.append(ch_name)

        self.missing_channels = []
        c = 0
        for ch_name in ch_names:
            for ptcp in self.participants:
                #raw = ptcp.raw_cleaned
                if ch_name not in ptcp.epochs_ransac_autoreject.info['ch_names']: #raw.info['ch_names']:
                    #print("Chan # ", c, "of ", len(ch_names), " missing")
                    self.missing_channels.append(ch_name)
            c+=1
        self.missing_channels = list(set(self.missing_channels))
        print("Missing channels ", self.missing_channels )

    def remove_missing_channels(self):
        r = 0
        for ptcp in self.participants:
            #raw = ptcp.raw_cleaned
            print("Remove missing channel loop: ptcp # ", r, " ID : ", ptcp.p_id)
            for missing_channel in self.missing_channels:
                status = missing_channel in ptcp.epochs_ransac_autoreject.info['ch_names']
                if status:
                    print("Dropping channel, ", missing_channel)
                    ptcp.epochs_ransac_autoreject.drop_channels([missing_channel])
                else:
                    print(missing_channel + " was not found, not removing")
            ptcp.signals_of_interest()
            r+=1
            
    def remove_missing_condition(self):
        ''' Remove participants who have the wrong condition for group analysis purposes
        '''
        
        self.identify_missing_channels()
        r = 0
        indices_to_delete = []
        for ptcp in self.participants:
            if conds_to_compare[events_to_tag[1]][0] not in ptcp.evoked_generic.keys() or conds_to_compare[events_to_tag[1]][1] not in ptcp.evoked_generic.keys():
                indices_to_delete.append(r)
                print("Remove missing condition loop : ptcp # ", r, " ID : ", ptcp.p_id, "Condition A: ", ptcp.cond_A, "Condition B: ", ptcp.cond_B, conds_to_compare[events_to_tag[1]][0], conds_to_compare[events_to_tag[1]][1])
            r+=1
        indices_to_delete = list(set(indices_to_delete))
        indices_to_delete = sorted(indices_to_delete, reverse=True)
        self.removed['Missing Condition'] = indices_to_delete
        print("Will delete based on missing condition ", indices_to_delete)
        for index in indices_to_delete:
            del self.participants[index]

    def remove_wrong_timespan(self):
        ''' Delete participants with an error due to inappropriate number of epochs
        '''
        r = 0
        indices_to_delete = []
        ptcps_to_delete = []
        for ptcp in self.participants:
            for condition in ptcp.evoked_generic.keys():
                L = len(ptcp.evoked_generic[condition].to_data_frame())
                #if L != NUM_EPOCH_SAMPLES:
                indices_to_delete.append(r)
                ptcps_to_delete.append(ptcp.p_id) 
                print("Remove wrong timespan loop: ptcp #", r , "Apparent Length: ", L, " Expected Length: ", NUM_EPOCH_SAMPLES)
            r+=1
        ptcps_to_delete = list(set(ptcps_to_delete))
        indices_to_delete = list(set(indices_to_delete))
        indices_to_delete = sorted(indices_to_delete, reverse=True)
        self.removed['Wrong Timespan'] = indices_to_delete        
        print("Will delete indices %s ptcps %s"%(str(indices_to_delete),str(ptcps_to_delete)))   
        for index in indices_to_delete:
            del self.participants[index]

    def remove_surplus_channels(self):
        '''Search and delete channel from any relevant object'''

        for ptcp in self.participants:
            required_channels = []
            for ch in range(1,NUM_CHANS_TO_KEEP+1+1): # MEG 001 to MEG 125 (MEG 071 will be missing in both child and adult systems, so 124 channels will be left over after this)
                ch_name = "MEG "+str(ch).zfill(3)
                required_channels.append(ch_name)
            ptcp.surplus_channels = []
            ptcp_chs = ptcp.epochs_ransac_autoreject.info['ch_names']
            for ptcp_ch in ptcp_chs:
                if ptcp_ch not in required_channels and "MEG" in ptcp_ch:
                    #print("Ptcp ID %s chan %s is surplus, dropping"%(str(ptcp.p_id), str(ptcp_ch)))
                    ptcp.surplus_channels.append(ptcp_ch)

            for attrib in ['raw_cleaned','evoked_generic','evoked_all','epochs','epochs_ransac','epochs_ransac_autoreject','epochs_ransac_autoreject_unstandardised']:
                if hasattr(ptcp,attrib):
                    obj = getattr(ptcp, attrib)
                    if type(obj) == type({}):
                        keys = list(obj.keys())
                        for key in keys:
                            for chan in ptcp.surplus_channels:
                                status = chan in obj[key].info['ch_names']
                                if status and "MEG" in chan:
                                    #print("Deleting channel 1 obj key, ", chan)
                                    obj[key].drop_channels([chan]) 
                                else:
                                    pass
                                    #print(chan + " was not found, not removing")
                    else:
                        for chan in ptcp.surplus_channels:
                            status = chan in obj.info['ch_names']
                            if status  and "MEG" in chan:
                                #print("Deleting channel 2 obj, ", chan)
                                obj.drop_channels([chan]) 
                            else:
                                pass
                                #print(chan + " was not found, not removing")            
     
    def rms_diff(self,condition='B'):
        self.squares_diff_grp = []
        num_chans = len(self.group_A_avg[condition].data)
        squares_A = []
        squares_B = []
        for time in range(0,len(self.group_A_avg[condition].data[0])): # Times
            sq_grpA_ttl = 0
            sq_grpB_ttl = 0  
            for chan in range(0,num_chans):
                sq_grpA = self.group_A_avg[condition].data[chan][time]**2 
                sq_grpB = self.group_B_avg[condition].data[chan][time]**2              
                    
                sq_grpA_ttl+=sq_grpA
                sq_grpB_ttl+=sq_grpB
            squares_A.append(sq_grpA)
            squares_B.append(sq_grpB)

            mean_sq_A = sq_grpA_ttl/num_chans    
            rms_A = math.sqrt(mean_sq_A)
            mean_sq_B = sq_grpB_ttl/num_chans    
            rms_B = math.sqrt(mean_sq_B) 
            sq_diff=rms_B-rms_A
            self.squares_diff_grp.append(sq_diff)
        print("Squares diff", self.squares_diff_grp)
        return self.squares_diff_grp

    def plot_rms_diff(self,condition):
        '''
        Plots the difference in GFPs between two groups and shows times where this is significant
        '''
        sample_ptcp = self.participants[0]
        times = sample_ptcp.epochs_ransac_autoreject[sample_ptcp.cond_A].times

        
        # Method A
        T_obs, clusters, cluster_p_values, H0 = mne.stats.permutation_cluster_1samp_test(
            np.array(self.squares_diff_grp),
            n_jobs=NUM_CPUS_Other,
            n_permutations=N_PERMUTATIONS,
            tail = 0,
            out_type='mask'
            #verbose='WARNING'
            #threshold=threshold,
            #adjacency=None,    
            )
        print("Clusters ", clusters)
        print("Cluster p values ", cluster_p_values)
        print("Mean ", np.mean(self.squares_diff_grp), "StD", np.std(self.squares_diff_grp))

        plt.plot(times, self.squares_diff_grp)
        plt.xlim([PLOT_START*SEC_TO_MS,PLOT_END*SEC_TO_MS])
        #plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed    
        plt.ylabel("GFP difference, fT", fontsize=FONTSIZE_LABELS)
        plt.title("Group B GFP minus Group A GFP, condition " +str(condition), fontsize=FONTSIZE_TITLE)
        plt.legend(["Group B - Group A"], fontsize=FONTSIZE_LEGEND, loc='best', bbox_to_anchor=(1, 0.5))

        # Put the cluster data in a viewable format
        width = 1
        p_clust = np.ones((width, width))
        for i_c, c in enumerate(clusters):
            c = c[0]
            #print(cluster_p_values[i_c])
            if cluster_p_values[i_c] <= CLUSTER_CUTOFF:
                h = plt.axvspan(times[c.start], times[c.stop - 1],
                            color='r', alpha=0.3)
            else:
                plt.axvspan(times[c.start], times[c.stop - 1], color=(0.3, 0.3, 0.3),
                            alpha=0.3)
                h = "No significant clusters"
        if len(clusters) == 0:
            h = [0,0.1]
        try:
            plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
            plt.ylabel("t values", fontsize=FONTSIZE_LABELS)
            #plt.legend(["all selfs cond B - cond A"], bbox_to_anchor=(0.75, 1.15), ncol=1, fontsize=FONTSIZE_LEGEND)
            plt.legend((h, ), ('cluster\n'+'p-value < '+str(CLUSTER_CUTOFF), ))
            hf = plt.plot(times, T_obs, 'g')    
        except Exception as e:
            pass
            #print("ERROR 3", e)
        plt.show()        


    def export_group(self,size_age_bucket_percent):

        '''Create a dataframe summarising differences between two groups
        '''
        df_gfps_group = pd.DataFrame({'time':           self.grand_average_cond_A_df['time'].values,                        
                                'GFP Group A EveryTrial':     self.rmses_A_everytrial,
                                'GFP Group B EveryTrial':     self.rmses_B_everytrial,
                                'GFP Group A EveryTrial minus Group B EveryTrial':   self.group_diff_everytrial,      # Expected to be +ve

                                'GFP Group A Cond A':   self.rmses_A_condA,
                                'GFP Group B Cond A':   self.rmses_B_condA,
                                'GFP Group A Cond A - Group B CondA':  self.group_diff_condA,       # Expected to be +ve

                                'GFP Group A Cond B':   self.rmses_A_condB,
                                'GFP Group B Cond B':   self.rmses_B_condB,
                                'GFP Group A Cond B - Group B CondB':  self.group_diff_condB,       # Expected to be -ve

                                'GFP Group B (B-A) minus Group A (B-A)': self.squares_diff_grp_mmr  # Expected to be -ve
                                      }
                                      )
        return df_gfps_group


    def group_analysis(self,
                        size_age_bucket_percent=PERCENTILES_AGE[1], # If e.g. 20, then there are 100/20 = 5 age groupings
                        age_bounds_low=[None,None],age_bounds_high=[None,None] # Specify e.g. [3,10], [15,20] if wanting to analyse those aged 3-10 vs those aged 15-20
                      ):
        '''
        Creates dataframes and images illustrating ERF and GFP differences between different age groups
        Saves these as .csvs and .jpegs.
        '''

        ## ERFs for each condition and MMR
        # Method A of getting age groupings
        age_num_groupings = int(100/size_age_bucket_percent)
        self.find_age_groupings(num_age_groupings=age_num_groupings)
        self.age_cutoff_low = self.age_bounds[0][1] # Code specific for num_age_groupings==2
        self.age_cutoff_high = self.age_bounds[age_num_groupings-1][0] # Code specific for num_age_groupings==2
        group_A_strings = []
        group_B_strings = []
        if age_bounds_low==[None,None] and age_bounds_high==[None,None]:
            for ptcp in self.participants:
                if ptcp.age <= self.age_cutoff_low:
                    group_A_strings.append(ptcp.p_id)
                if ptcp.age > self.age_cutoff_high:
                    group_B_strings.append(ptcp.p_id)
        # Method B of getting age groupings            
        else:
            assert(age_bounds_low[1]) < age_bounds_high[0], "Choose valid age ranges, categories cannot share members"
            assert(age_bounds_low[0]) < age_bounds_low[1], "Choose lower bound of age that's lower than high bound"
            assert(age_bounds_high[0]) < age_bounds_high[1], "Choose lower bound of age that's lower than high bound"        
            for ptcp in self.participants:
                if ptcp.age >=age_bounds_low[0] and ptcp.age <= age_bounds_low[1]:
                    group_A_strings.append(ptcp.p_id)
                if ptcp.age >=age_bounds_high[0] and ptcp.age <= age_bounds_high[1]:
                    group_B_strings.append(ptcp.p_id)                
                    
        self.group_A, self.group_B = group_A_strings, group_B_strings
        print("Group A: ", self.group_A, "Group B: ", self.group_B)
        
        # Group 1 grand average
        self.group_A_avg = {}
        self.group_B_avg = {}

        self.group_A_avg[None] = self.average_evoked(self.group_A,condition=None) # Get the average for condition B amongst every member of Group A 
        self.group_A_avg['B']  = self.average_evoked(self.group_A,condition='B') # Get the average for condition B amongst every member of Group A        
        self.group_A_avg['A']  = self.average_evoked(self.group_A,condition='A') # Get the average for condition A amongst every member of Group A
        self.group_B_avg[None] = self.average_evoked(self.group_B,condition=None) # Get the average for condition B amongst every member of Group A     
        self.group_B_avg['B']  = self.average_evoked(self.group_B,condition='B') # Get the average for condition B amongst every member of Group B
        self.group_B_avg['A']  = self.average_evoked(self.group_B,condition='A') # Get the average for condition A amongst every member of Group B
       # Convert to DFs
        self.group_A_avg_condB_df = self.group_A_avg['B'].to_data_frame()
        self.group_A_avg_condA_df = self.group_A_avg['A'].to_data_frame()
        self.group_B_avg_condB_df = self.group_B_avg['B'].to_data_frame()
        self.group_B_avg_condA_df = self.group_B_avg['A'].to_data_frame()

        timestring = str(int(self.participants[0].epochs_ransac_autoreject.times[0]*SEC_TO_MS+EPSILON_TIME))+ "-"+str(int(self.participants[0].epochs_ransac_autoreject.times[-1]*SEC_TO_MS+EPSILON_TIME))+"ms"
        filename = "ERF Group A Cond EveryTrial "+str(size_age_bucket_percent)+"%ile"+" "+timestring+".csv"
        self.group_A_avg[None].to_data_frame().to_csv(filename,index=False)       
        filename = "ERF Group A Cond A "+str(size_age_bucket_percent)+"%ile"+" "+timestring+".csv"
        self.group_A_avg['A'].to_data_frame().to_csv(filename,index=False)        
        filename = "ERF Group A Cond B "+str(size_age_bucket_percent)+"%ile"+" "+timestring+".csv"
        self.group_A_avg['B'].to_data_frame().to_csv(filename,index=False)    
        filename = "ERF Group B Cond EveryTrial "+str(size_age_bucket_percent)+"%ile"+" "+timestring+".csv"
        self.group_B_avg[None].to_data_frame().to_csv(filename,index=False)       
        filename = "ERF Group B Cond A "+str(size_age_bucket_percent)+"%ile"+" "+timestring+".csv"
        self.group_B_avg['A'].to_data_frame().to_csv(filename,index=False)        
        filename = "ERF Group B Cond B "+str(size_age_bucket_percent)+"%ile"+" "+timestring+".csv"
        self.group_B_avg['B'].to_data_frame().to_csv(filename,index=False)    


        self.group_A_avg_MMR_df = self.group_A_avg_condB_df - self.group_A_avg_condA_df
        self.group_B_avg_MMR_df = self.group_B_avg_condB_df - self.group_B_avg_condA_df 
        # Somewhat meaningless 'sum' column       
        self.group_A_avg_MMR_df = sum_df(self.group_A_avg_MMR_df,mode='plain')  
        self.group_B_avg_MMR_df = sum_df(self.group_B_avg_MMR_df,mode='plain')              
        # Plot MMRs
        self.group_A_avg_MMR_df['time'] = self.group_A_avg_condB_df['time'].values
        self.group_B_avg_MMR_df['time'] = self.group_B_avg_condB_df['time'].values
        
        ## Calculate difference in the MMR between groups
        # Method A
        self.group_diff_MMR_df = self.group_A_avg_MMR_df - self.group_B_avg_MMR_df
        self.group_diff_MMR_df['time'] = self.group_A_avg_MMR_df['time']

        # Plot difference in the MMR between groups
        # plt.plot(self.group_diff_MMR_df['time'].values, self.group_diff_MMR_df['sum'].values)
        # plt.show()
        filename = "ERF MMR Group A minus MMR Group B "+str(size_age_bucket_percent)+"%"+" "+timestring+".csv"
        self.group_diff_MMR_df.to_csv(filename,index=False)       


        # GFP condition EVERY TRIAL for both groups
        self.squares_diff_grp_cond_everytrial = []
        num_chans = len(self.group_A_avg[None].data)
        squares_A = []
        squares_B = []
        self.rmses_A_everytrial = []
        self.rmses_B_everytrial = []   
        for time in range(0,len(self.group_A_avg[None].data[0])): # Times
            sq_grpA_ttl = 0
            sq_grpB_ttl = 0  
            for chan in range(0,num_chans):
                sq_grpA = self.group_A_avg[None].data[chan][time]**2 # self.average_evoked(self.group_A)[chan][time]**2
                sq_grpA_ttl+=sq_grpA
                sq_grpB = self.group_B_avg[None].data[chan][time]**2 
                sq_grpB_ttl+=sq_grpB
            squares_A.append(sq_grpA)
            squares_B.append(sq_grpB)
            mean_sq_A = sq_grpA_ttl/num_chans    
            rms_A = math.sqrt(mean_sq_A)
            mean_sq_B = sq_grpB_ttl/num_chans    
            rms_B = math.sqrt(mean_sq_B) 
            self.rmses_A_everytrial.append(rms_A)
            self.rmses_B_everytrial.append(rms_B)        
            rms_diff=rms_B-rms_A
            self.squares_diff_grp_cond_everytrial.append(rms_diff)
            
        self.group_A_everytrial             = RMS_DF(self.group_A_avg[None].to_data_frame())
        self.group_B_everytrial             = RMS_DF(self.group_B_avg[None].to_data_frame())
        self.group_diff_everytrial          = [(self.group_A_everytrial[x] - self.group_B_everytrial[x]) for x in range(0,len(self.group_A_everytrial))]
        self.group_diff_everytrial_df       = pd.DataFrame({'time': self.grand_average_cond_A_df['time'].values, 'GFP_DIFF': self.group_diff_everytrial})

        filename = "GFP Group B every trial minus Group A every trial "+str(size_age_bucket_percent)+"%"+" "+timestring+".csv" # Names deliberately swapped
        self.group_diff_everytrial_df.to_csv(filename,index=False)   


        # GFP condition A for both groups
        self.squares_diff_grp_cond_A = []
        num_chans = len(self.group_A_avg['A'].data)
        squares_A = []
        squares_B = []
        self.rmses_A_condA = []
        self.rmses_B_condA = []       
        for time in range(0,len(self.group_A_avg['A'].data[0])): # Times
            sq_grpA_ttl = 0
            sq_grpB_ttl = 0  
            for chan in range(0,num_chans):
                sq_grpA = self.group_A_avg['A'].data[chan][time]**2 # self.average_evoked(self.group_A)[chan][time]**2
                sq_grpA_ttl+=sq_grpA
                sq_grpB = self.group_B_avg['A'].data[chan][time]**2 
                sq_grpB_ttl+=sq_grpB
            squares_A.append(sq_grpA)
            squares_B.append(sq_grpB)
            mean_sq_A = sq_grpA_ttl/num_chans    
            rms_A = math.sqrt(mean_sq_A)
            mean_sq_B = sq_grpB_ttl/num_chans    
            rms_B = math.sqrt(mean_sq_B) 
            self.rmses_A_condA.append(rms_A)
            self.rmses_B_condA.append(rms_B)          
            rms_diff=rms_B-rms_A
            self.squares_diff_grp_cond_A.append(rms_diff)

        self.group_A_condA            = RMS_DF(self.group_A_avg['A'].to_data_frame())
        self.group_B_condA            = RMS_DF(self.group_B_avg['A'].to_data_frame())
        self.group_diff_condA         = [self.group_A_condA[x] - self.group_B_condA[x] for x in range(0,len(self.group_A_condA))]
        self.group_diff_condA_df      = pd.DataFrame({'time': self.grand_average_cond_A_df['time'].values, 'GFP_DIFF_COND_A': self.group_diff_condA})    
        
    
        # GFP condition B for both groups
        self.squares_diff_grp_cond_B = []
        num_chans = len(self.group_A_avg['B'].data)
        squares_A = []
        squares_B = []
        self.rmses_A_condB = []
        self.rmses_B_condB = []    
        for time in range(0,len(self.group_A_avg['B'].data[0])): # Times
            sq_grpA_ttl = 0
            sq_grpB_ttl = 0  
            for chan in range(0,num_chans):
                sq_grpA = self.group_A_avg['B'].data[chan][time]**2 # self.average_evoked(self.group_A)[chan][time]**2
                sq_grpA_ttl+=sq_grpA
                sq_grpB = self.group_B_avg['B'].data[chan][time]**2 # self.average_evoked(self.group_B)[chan][time]**2
                sq_grpB_ttl+=sq_grpB
            squares_A.append(sq_grpA)
            squares_B.append(sq_grpB)
                
            mean_sq_A = sq_grpA_ttl/num_chans    
            rms_A = math.sqrt(mean_sq_A)
            mean_sq_B = sq_grpB_ttl/num_chans    
            rms_B = math.sqrt(mean_sq_B) 
            self.rmses_A_condB.append(rms_A)
            self.rmses_B_condB.append(rms_B)          
            rms_diff=rms_A-rms_B        # Old minus young
            self.squares_diff_grp_cond_B.append(rms_diff)

        self.group_A_condB            = RMS_DF(self.group_A_avg['B'].to_data_frame())
        self.group_B_condB            = RMS_DF(self.group_B_avg['B'].to_data_frame())
        self.group_diff_condB         = [self.group_A_condB[x] - self.group_B_condB[x] for x in range(0,len(self.group_A_condB))]
        self.group_diff_condB_df      = pd.DataFrame({'time': self.grand_average_cond_A_df['time'].values, 'GFP_DIFF_COND_B': self.group_diff_condB})    
        
        # GFP "MMR" (condition B GFP minus condition A GFP) for group B minus same for group A
        self.squares_diff_grp_mmr     = [self.squares_diff_grp_cond_B[x] - self.squares_diff_grp_cond_A[x]  for x in range(0,len(self.squares_diff_grp_cond_A))]
        self.group_diff_MMR_GFP_df    = pd.DataFrame({'time': self.grand_average_cond_A_df['time'].values, 'GFP_DIFF_MMR': self.squares_diff_grp_mmr})
        

    def save_gfp_diff_cond(self,file_name='Normal',filter_age_max = None):
        '''
        Tests difference between GFPs based on condition
        Uses individuals' average evoked response each as inputs to the cluster analysis (which is 74 participants x 100 time-points)
        '''
        all_gfps_cond_A = []
        all_gfps_cond_B = []
        for ptcp in self.participants:
            p_id = str(ptcp.p_id)
            avg_condition_A = ptcp.epochs_ransac_autoreject[ptcp.cond_A].average().to_data_frame()
            avg_condition_B = ptcp.epochs_ransac_autoreject[ptcp.cond_B].average().to_data_frame()
            gfp_condition_A = np.array(RMS_DF(avg_condition_A))
            gfp_condition_B = np.array(RMS_DF(avg_condition_B))
            

            if filter_age_max!=None:
                extra_str = str(filter_age_max)
                title = 'GFP in high surprise minus\n'+'GFP in low surprise trials\n'+"Max age "+str(extra_str)
                filename = "GFP vs condition, differences "+str(int(100/NUM_BINS_SURPRISE))+"%, "+"Max age "+str(extra_str)+" Sig test" if file_name=='Normal' else None
                
                if ptcp.age < filter_age_max:
                    all_gfps_cond_A.append(gfp_condition_A)
                    all_gfps_cond_B.append(gfp_condition_B)
            else:
                extra_str = ""
                title = 'GFP in high surprise minus\n'+'GFP in low surprise trials\n'
                filename = "GFP vs condition, differences "+str(int(100/NUM_BINS_SURPRISE))+"%, "+" Sig test" if file_name=='Normal' else None
                all_gfps_cond_A.append(gfp_condition_A)
                all_gfps_cond_B.append(gfp_condition_B)
                
                    
        all_gfps_cond_A = np.array(all_gfps_cond_A)
        all_gfps_cond_B = np.array(all_gfps_cond_B)



        #print("TIMES HERE ARE :", self.grand_average_cond_A_df['time'].values)
        self.cluster_analysis_new(          conditionA = all_gfps_cond_A,
                                            conditionB = all_gfps_cond_B,
                                            compare_mode = "Condition",
                                            filename = filename,
                                            title =  title, # 'GFP Cond B - GFP Cond A',
                                            legend = ["GFP Cond B - Cond A"],
                                            times =  self.grand_average_cond_A_df['time'].values,
                                            gfp_mode = True, group_mode=True, extract_data=False,save_cohens_d=True) 
    
    def save_gfp_diff_ages(self, size_age_bucket_percent=10, override_filename=None):
        '''
        Tests difference between GFPs based on condition X age grouping
        Uses individuals' average evoked response each as inputs to the cluster analysis (which is 74 participants x 100 time-points)
        '''    
        
        all_gfps_grp_A_cond_everytrial = []
        all_gfps_grp_B_cond_everytrial = []    
        all_gfps_grp_A_cond_A = []
        all_gfps_grp_B_cond_A = []
        all_gfps_grp_A_cond_B = []
        all_gfps_grp_B_cond_B = []
        # The 'MMR' of the GFP per age group is (GFP cond B - GFP cond A)
        all_gfps_grp_A_MMR    = []
        all_gfps_grp_B_MMR    = []
        
        for ptcp in self.participants:
            p_id = str(ptcp.p_id)
            avg_condition_everytrial = ptcp.epochs_ransac_autoreject.average().to_data_frame()        
            avg_condition_A = ptcp.epochs_ransac_autoreject[ptcp.cond_A].average().to_data_frame()
            avg_condition_B = ptcp.epochs_ransac_autoreject[ptcp.cond_B].average().to_data_frame()
            gfp_condition_everytrial = np.array(RMS_DF(avg_condition_everytrial))  
            
            rms_cond_A = RMS_DF(avg_condition_A)  # Gfp of the grand average for condition BA
            gfp_condition_A = np.array(rms_cond_A)
            rms_cond_B = RMS_DF(avg_condition_B) # Gfp of the grand average for condition B
            gfp_condition_B = np.array(rms_cond_B)
            gfp_mmr = np.array( [rms_cond_B[x] - rms_cond_A[x] for x in range(0,len(rms_cond_A)) ] ) # Difference in GFP
                  
            if ptcp.p_id in self.group_A:
                all_gfps_grp_A_cond_everytrial.append(gfp_condition_everytrial)            
                all_gfps_grp_A_cond_A.append(gfp_condition_A)
                all_gfps_grp_A_cond_B.append(gfp_condition_B)
                all_gfps_grp_A_MMR.append(gfp_mmr) 
            if ptcp.p_id in self.group_B:
                all_gfps_grp_B_cond_everytrial.append(gfp_condition_everytrial)            
                all_gfps_grp_B_cond_A.append(gfp_condition_A)
                all_gfps_grp_B_cond_B.append(gfp_condition_B)    
                all_gfps_grp_B_MMR.append(gfp_mmr)
                
        all_gfps_grp_A_cond_everytrial = np.array(all_gfps_grp_A_cond_everytrial)
        all_gfps_grp_B_cond_everytrial = np.array(all_gfps_grp_B_cond_everytrial)
        all_gfps_grp_A_cond_A = np.array(all_gfps_grp_A_cond_A)
        all_gfps_grp_B_cond_A = np.array(all_gfps_grp_B_cond_A)
        all_gfps_grp_A_cond_B = np.array(all_gfps_grp_A_cond_B)
        all_gfps_grp_B_cond_B = np.array(all_gfps_grp_B_cond_B)
        all_gfps_grp_A_MMR = np.array(all_gfps_grp_A_MMR)
        all_gfps_grp_B_MMR = np.array(all_gfps_grp_B_MMR)    
        
        if override_filename!=None:
            size_age_bucket_percent = override_filename
        else:
            size_age_bucket_percent = str(size_age_bucket_percent)+"%"+"ile"

        timestring = str(int(self.participants[0].epochs_ransac_autoreject.times[0]*SEC_TO_MS+EPSILON_TIME))+ "-"+str(int(self.participants[0].epochs_ransac_autoreject.times[-1]*SEC_TO_MS+EPSILON_TIME))+"ms"
        self.cluster_analysis_new(          conditionA = all_gfps_grp_B_cond_everytrial, # This gets switched around, so we subtract GFP for Group A from Group B (diff should be > 0)
                                            conditionB = all_gfps_grp_A_cond_everytrial,
                                            compare_mode = "Age",                                            
                                            filename = 'GFP vs age (every trial, '+size_age_bucket_percent+'), sig test',
                                            title = 'GFP, every trial\n'+'Young minus old ('+size_age_bucket_percent+' age split)',
                                            legend = ["All conditions\n"+"Young minus old"],
                                            times = self.grand_average_cond_A_df['time'].values,
                                            gfp_mode = True, group_mode=True, extract_data=False) 
        
        self.cluster_analysis_new(          conditionA = all_gfps_grp_B_cond_A,
                                            conditionB = all_gfps_grp_A_cond_A,
                                            compare_mode = "Age",
                                            filename = 'GFP vs age '+size_age_bucket_percent+' (Low surprise, '+str(100/NUM_BINS_SURPRISE)+'%), sig test',
                                            title = 'GFP, low surprise\n'+'Young minus old  ('+size_age_bucket_percent+' age split)',
                                            legend = ["Low surprise\n"+"Young minus old"],
                                            times = self.grand_average_cond_A_df['time'].values,
                                            gfp_mode = True, group_mode=True, extract_data=False)     

        self.cluster_analysis_new(          conditionA = all_gfps_grp_B_cond_B,
                                            conditionB = all_gfps_grp_A_cond_B,
                                            compare_mode = "Age",
                                            filename = 'GFP vs age '+size_age_bucket_percent+' (High surprise, '+str(100/NUM_BINS_SURPRISE)+'%), sig test',
                                            title = 'GFP, high surprise\n'+'Young minus old ('+size_age_bucket_percent+' age split)',
                                            legend = ["High surprise\n"+"Young minus old"],
                                            times = self.grand_average_cond_A_df['time'].values,
                                            gfp_mode = True, group_mode=True, extract_data=False)  
        
        self.cluster_analysis_new(          conditionA = all_gfps_grp_B_MMR,
                                            conditionB = all_gfps_grp_A_MMR,
                                            compare_mode = "Age",
                                            filename = 'GFP vs age (MMR, '+str(size_age_bucket_percent)+'), sig test',
                                            title = '(GFP high surprise - GFP low surprise)'+', '+str(100/NUM_BINS_SURPRISE)+'% cond split\n'+'Young minus old ('+str(size_age_bucket_percent)+' age split)',
                                            legend = ["aMMR\n"+"Young minus old"],
                                            times = self.grand_average_cond_A_df['time'].values,
                                            gfp_mode = True, group_mode=True, extract_data=False)      


    def get_child_adult_sample(self):

        '''Find any child in the group for use of their headshape'''
        if not hasattr(self,'sample_child'):
            found_child = False
            for ptcp in self.participants:
                if not found_child:
                    if not ptcp.is_adult_system:
                        self.sample_child = copy.deepcopy(ptcp)
                        found_child = True
        if not hasattr(self,'sample_adult'):
            found_adult = False
            for ptcp in self.participants:
                if not found_adult:
                    if ptcp.is_adult_system:
                        self.sample_adult = copy.deepcopy(ptcp)
                        found_adult = True                        


    def cluster_analysis_new(self, conditionA, conditionB, times, filename=None, title=None, legend = None,
                                gfp_mode=False, group_mode=False, extract_data=True,
                                condition_label = "B", n_permutations=N_PERMUTATIONS, cluster_cutoff=CLUSTER_CUTOFF,
                                compare_mode = "Age", # 'Age' or 'Condition'
                                save_cohens_d = True):

        '''Perform cluster analysis to examine for significant differences in the ERF or GFP (depending on choice)

        Typically use
        -> Group B minus Group A
        -> Condition B minus Condition A
        Describe which with compare_mode
        
        Other params
            save_cohens_d: whether to save these to the .csv
            group_mode, condition_label: used for setting title (if title not set)
            extract_data: describes the format of the data input
            gfp_mode: If True, we are comparing GFPs. If False, we are comparing ERFs
        '''

        timepoints = self.participants[0].epochs_ransac_autoreject.times
        timestring = str(int(timepoints[0]*SEC_TO_MS+EPSILON_TIME))+ "-"+str(int(timepoints[-1]*SEC_TO_MS+EPSILON_TIME))+"ms"


        if extract_data:
            conditionA = conditionA.data
            conditionB = conditionB.data

        self.get_child_adult_sample()
        evok_avg = self.sample_child.epochs_ransac_autoreject[self.sample_child.cond_A].average()


        info = evok_avg.info  
        #headshape = copy.deepcopy(self.group_A_avg['B'])
        headshape = self.sample_child.epochs_ransac_autoreject[self.sample_child.cond_A]
        #headshape = centre_sensor_locations(headshape)[0]
        info = headshape.info
        # Adjacency matrix
        ch_adjacency, ch_names = mne.channels.find_ch_adjacency(info, 'mag');

        X = np.array([conditionB, conditionA])
        # print(np.shape(X))
        T_obs, clusters, cluster_p_values, H0 = mne.stats.permutation_cluster_test(X,   
                                                                                    n_permutations=n_permutations,      #permutation_cluster_test([conditionB, conditionA], n_permutations=n_permutations,
                                                                                    tail=0,
                                                                                    n_jobs=NUM_CPUS_Other,
                                                                                    out_type='mask',
                                                                                    verbose='DEBUG',
                                                                                    #threshold=threshold, 
                                                                                    #adjacency=ch_adjacency
                                                                                    )


        #Plot
        print("T obs: ", T_obs)
        print("Cluster p values ", cluster_p_values)

        from matplotlib import pyplot as plt
        try:
            plt.close('all')
            plt.clf()
        except:
            print("Problem closing graphs")
        plt.subplot(211)
        if group_mode:
            if gfp_mode:
                if title==None:
                    plt.title('GFP for Group B (older) minus group A (younger), condition='+str(condition_label), fontsize=FONTSIZE_TITLE) 
                else:
                    plt.title(title, fontsize=FONTSIZE_TITLE)
                plt.ylabel("Difference in GFP", fontsize=FONTSIZE_LABELS)
            else:
                if title==None:
                    plt.title('ERF for Group B (older) minus group A (younger), condition='+str(condition_label), fontsize=FONTSIZE_TITLE) 
                else:
                    plt.title(title, fontsize=FONTSIZE_TITLE)
                plt.ylabel("GFP difference", fontsize=FONTSIZE_LABELS)
        else:
            if gfp_mode:
                if title==None:            
                    plt.title("GFP for high surprise trials (condition A)\nminus GFP for low surprise trials, condition="+str(condition_label), fontsize=FONTSIZE_TITLE)
                else:
                    plt.title(title, fontsize=FONTSIZE_TITLE)
                plt.ylabel("Difference in GFP", fontsize=FONTSIZE_LABELS)
            else:
                if title==None:            
                    plt.title("ERF for high surprise trials (condition A)\nminus ERF for low surprise trials, condition="+str(condition_label), fontsize=FONTSIZE_TITLE)
                else:
                    plt.title(title, fontsize=FONTSIZE_TITLE)
                plt.ylabel("MEG signal difference (fT)", fontsize=FONTSIZE_LABELS)

        diff = conditionB.mean(axis=0) - conditionA.mean(axis=0)
        if len(diff) != len(times):
            self.warn("DIFFERENT LENGTH MATRICES TO NUMBER OF TIME SNAPSHOTS") # Then assumes our analysis begins at the start of the epoch
            diff = diff[0:len(times)]
            T_obs = T_obs[0:len(times)]

        plt.plot(times, diff)
        plt.ylabel("Difference\n"+"(fT)", fontsize=FONTSIZE_LABELS)        
        #plt.ylim([-12,6])
        plt.xlim([PLOT_START*SEC_TO_MS,PLOT_END*SEC_TO_MS])
        #plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed    

        if legend == None:
            if group_mode:
                plt.legend(["Group B - Group A"], fontsize=FONTSIZE_LEGEND, loc='best', bbox_to_anchor=(1, 0.5))            
            else:
                plt.legend(["Cond B - Cond A"], fontsize=FONTSIZE_LEGEND, loc='best', bbox_to_anchor=(1, 0.5))
        else:
            plt.legend(legend, fontsize=FONTSIZE_LEGEND, loc='best', bbox_to_anchor=(1, 0.5))
            
        plt.subplot(212)
        times = list(times)
        #print("TIMES ARE 2", times)
        for i_c, c in enumerate(clusters):
            c = c[0]
            if cluster_p_values[i_c] <= cluster_cutoff:
                h = plt.axvspan(times[c.start], times[c.stop - 1],
                                color='r', alpha=0.3)
            else:
                plt.axvspan(times[c.start], times[c.stop - 1], color=(0.3, 0.3, 0.3),
                            alpha=0.3)
                h = "No significant clusters"
        if len(clusters) == 0:
            h = [0,0.1]
        hf = plt.plot(times, T_obs, 'g')
        plt.xlabel("Time (ms)", fontsize = FONTSIZE_LABELS)
        plt.ylabel("t values", fontsize = FONTSIZE_LABELS)
        #plt.ylim([0,20])

        #plt.legend(["all Grps cond B - cond A"], bbox_to_anchor=(0.75, 1.15), ncol=1, fontsize=FONTSIZE_LEGEND)
        plt.legend((h, ), ('cluster\n'+'p-value < '+str(cluster_cutoff), ), fontsize = FONTSIZE_LEGEND)
        if filename!=None:
            plt.savefig(filename+ " "+timestring+".jpg",
                            format='jpeg',
                            dpi=DPI,
                            bbox_inches='tight')
            EPOCH_SIZE_MS = int((times[-1]-times[0])/   ((len(times)-1)   ) + EPSILON_TIME)
            with open(filename+" "+timestring+' p&t values'+".txt", 'w') as f:
                f.write("As slices: " + str(list(clusters)))
                f.write("\n")
                f.write("Times:\n") 
                for cluster in list(clusters):
                    cluster = cluster[0]
                    f.write(str(int(times[0])+  int(cluster.start*EPOCH_SIZE_MS))+"ms to "+str(int(times[0])+  int((cluster.stop-1)*EPOCH_SIZE_MS))+ " ms, ")
                f.write("\n")
                f.write("p values:\n")# +str(list(cluster_p_values))) 
                for i_c_inner, c_inner in enumerate(clusters): # Writes one p-value per cluster (shaded region)
                    f.write(str(cluster_p_values[i_c_inner])+", ") 
                f.write("\n")
                f.write("t values: " +str(T_obs))
                f.write("\n")
                f.write("\n")
                f.write("t values at sig times\n")
                for cluster in list(clusters):
                    cluster = cluster[0]
                    print("Start: ", cluster.start, " End: ", cluster.stop)
                    for n in range(cluster.start,cluster.stop):
                        print("n ", n)
                        f.write(str(T_obs[n])+", ")
                    f.write("\n")



            df_gfp_diffs                = pd.DataFrame()
            df_gfp_diffs['time']        = times
            col_names                   = ["Low surprise (A)", "High surprise (B)"] if compare_mode == 'Condition' else ["Old (A)","Young (B)"]
            col_name_diff               = "High (B) - Low (A)" if compare_mode == 'Condition' else "Young (B) - Old (A)"
            df_gfp_diffs[col_names[0]]  = conditionA.mean(axis=0)
            df_gfp_diffs[col_names[1]]  = conditionB.mean(axis=0)
            df_gfp_diffs[col_name_diff] = diff
            df_gfp_diffs['t-values']    = T_obs
            if save_cohens_d:
                c = 0 # Count for condition
                for gfps in [conditionA,conditionB]:
                    cond_code = 'A' if c==0 else 'B'
                    gfps_timepoint = {}
                    t = 0 # Counter when sweeping through time
                    for t in range(0,len(timepoints)):
                        p = 0 # Counter when sweeping through participants
                        for ptcp in gfps:
                            #print(p,t,key)
                            if timepoints[t] in gfps_timepoint.keys():
                                gfps_timepoint[timepoints[t]].append(gfps[p][t])
                            else:
                                gfps_timepoint[timepoints[t]] = [gfps[p][t]]
                            p+=1
                        t+=1

                    means = []
                    stds = []
                    for time in sorted(gfps_timepoint.keys()):
                        mean = np.mean(gfps_timepoint[time])
                        std = np.std(gfps_timepoint[time])
                        print("Time %s mean %s std %s"%(str(time),str(mean),str(std)))
                        means.append(mean)
                        stds.append(std)
                    df_gfp_diffs['mean_cond_'+str(cond_code)] = means
                    df_gfp_diffs['std_cond_'+str(cond_code)] = stds
                    c+=1

                denominators = [    math.sqrt((df_gfp_diffs['std_cond_A'].values[x]**2+df_gfp_diffs['std_cond_B'].values[x]**2)/2) for x in range(0,len(df_gfp_diffs['std_cond_B'].values))    ]
                cohens_ds_1    = [ df_gfp_diffs[col_name_diff].values[x]/denominators[x] for x in range(0,len(df_gfp_diffs[col_name_diff].values)) ] 
                cohens_ds_2    = [ (df_gfp_diffs['mean_cond_B'].values[x] - df_gfp_diffs['mean_cond_A'].values[x])/denominators[x] for x in range(0,len(df_gfp_diffs['mean_cond_A'].values)) ] 
                                
                df_gfp_diffs['cohens_d_1'] = cohens_ds_1
                df_gfp_diffs['cohens_d_2'] = cohens_ds_1
            df_gfp_diffs.to_csv(filename+ ' ' + timestring+".csv",index=False) 

        plt.show()
        return diff     

    def return_method_dct(self,str_):
        '''Return a dictionary with a key's value modified to True
        '''
        method_dct = {  'On every trial':                        False,     # Do regression over all trials
                        'On high surprise trials':              False,     # Only do regression over 'deviant'/surprising/condition_B *events*.                                             
                        'On low surprise trials':               False,     # Only do regression over 'pre-deviant'/unsurprising/condition_A *events*. 

                        'Subtracting average of all trials':    False,     # Subtract average from *every* epoch 
                        'Subtracting high surprise average':    False,     # Subtract high surprise average from every high surprise epoch
                        'Subtracting low surprise average':     False,     # Subtract low surprise average from every low surprise epoch
                        'Subtracting deviants only':            False      # Subtract 'deviants'/surprising/condition_B from 'predeviants'/less-surprising/condition_A
                     }

        method_dct[str_] = True
        return method_dct
        
    def return_waveform_str(self,method_dct):
        '''Produce a string to describe the waveform used to calculate the data upon which to calculate a statistic (E.g. correlation) between predictor and ERF/GFP
        '''
        for key in method_dct.keys():
            if method_dct[key] == True:
                return key
        return ""      

    def cluster_analysis_spatio_temporal_v2(self,mode='condition',condition_comparison='B', size_age_bucket_percent=20):
        '''
        Produces graphs of where the spatiotemporal clusters exist for ERFs
        a) Differences between conditions
        b) Differences between age groups
        c) MMR differences per age
        d) Age by condition interactions
        '''

        self.remove_surplus_channels()
        self.identify_missing_channels()
        self.get_child_adult_sample()
        condition_A = self.sample_child.cond_A
        condition_B = self.sample_child.cond_B
        condition_names = [condition_A,condition_B]    
        self.compare_ERFs_conditions(cond_A=condition_A,cond_B=condition_B);
        REF_OBJECT = mne.combine_evoked([self.group_B_avg['B'],self.group_B_avg['A'], self.group_A_avg['B'],self.group_A_avg['A']], weights=[-1,1,1,-1])
        chans_data = REF_OBJECT.data
        times_data = REF_OBJECT.data[0]
 
        # Evokeds
        evoked = {}
        all_evoked_condition_A = []
        all_evoked_condition_B = [] 
        for ptcp in self.participants:
            evoked[ptcp.p_id] = {}        
            ptcp.epochs_ransac_autoreject = ptcp.epochs_ransac_autoreject.reorder_channels(sorted(ptcp.epochs_ransac_autoreject.ch_names))
        
        if mode == 'condition':
            for ptcp in self.participants:
                p_id = str(ptcp.p_id)
                avg_condition_A = ptcp.epochs_ransac_autoreject[ptcp.cond_A].average()  # ptcp.evoked_generic[ptcp.cond_A] # 
                avg_condition_B = ptcp.epochs_ransac_autoreject[ptcp.cond_B].average()  # ptcp.evoked_generic[ptcp.cond_B] # 
                evoked[p_id][ptcp.cond_A] = avg_condition_A
                all_evoked_condition_A.append(avg_condition_A)
                evoked[p_id][ptcp.cond_B] = avg_condition_B
                all_evoked_condition_B.append(avg_condition_B)
        else:
            age_num_groupings = int(100/size_age_bucket_percent)
            self.find_age_groupings(num_age_groupings=age_num_groupings)
            age_cutoff_low = self.age_bounds[0][1] # Code specific for num_age_groupings==2
            age_cutoff_high = self.age_bounds[age_num_groupings-1][0] # Code specific for num_age_groupings==2
            print("Age number of groupings %s cutoff low %s cutoff high %s"%(str(age_num_groupings),str(age_cutoff_low),str(age_cutoff_high)))

            if mode == "age": # Compare two age groupings over all conditions 
                for ptcp in self.participants:
                    if ptcp.age <= age_cutoff_low:
                        all_evoked_condition_A.append(ptcp.epochs_ransac_autoreject.average())  #  Group A(Youngest x% percentile) - if standardisation occurs, this is done per condition
                    elif ptcp.age >= age_cutoff_high:
                        all_evoked_condition_B.append(ptcp.epochs_ransac_autoreject.average())  #  Group B(Oldest x% percentile)
            if mode == "mmr_by_age":
                for ptcp in self.participants:
                    if ptcp.age <= age_cutoff_low:
                        mmr = mne.combine_evoked([ptcp.evoked_generic[ptcp.cond_B],ptcp.evoked_generic[ptcp.cond_A]],[1,-1])
                        all_evoked_condition_A.append(mmr) # Group A(Youngest x% percentile) condition A
                    elif ptcp.age >= age_cutoff_high:
                        mmr = mne.combine_evoked([ptcp.evoked_generic[ptcp.cond_B],ptcp.evoked_generic[ptcp.cond_A]],[1,-1])
                        all_evoked_condition_B.append(mmr) # Group A(Youngest x% percentile) condition A
            if mode == "condition_by_age": #  # Compare two age groupings over each condition separately
                for ptcp in self.participants:
                    if ptcp.age <= age_cutoff_low:
                        if condition_comparison == "A":
                            all_evoked_condition_A.append(ptcp.epochs_ransac_autoreject[ptcp.cond_A].average()) # Group A(Youngest x% percentile) condition A
                        elif condition_comparison == "B":
                            all_evoked_condition_A.append(ptcp.epochs_ransac_autoreject[ptcp.cond_B].average()) # Group A(Youngest x% percentile) condition B
                    elif ptcp.age >= age_cutoff_high:
                        if condition_comparison == "A":
                            all_evoked_condition_B.append(ptcp.epochs_ransac_autoreject[ptcp.cond_A].average()) # Group B(Oldest x% percentile) condition A
                        elif condition_comparison == "B":
                            all_evoked_condition_B.append(ptcp.epochs_ransac_autoreject[ptcp.cond_B].average()) # Group B(Oldest x% percentile) condition B


        #evokeds = [mne.combine_evoked([evoked_A for evoked_A in all_evoked_condition_A],weights=[1 for evoked_A in all_evoked_condition_A]), mne.combine_evoked([evoked_B for evoked_B in all_evoked_condition_B],weights=[1 for evoked_B in all_evoked_condition_B])];
        min_dim = min(len(all_evoked_condition_A), len(all_evoked_condition_B))
        multiplier_cond_A, multiplier_cond_B = 1, 1
        if len(all_evoked_condition_A)!=len(all_evoked_condition_B):
            # Require same dimension
            print("Mismatched lengths cond A %s cond B %s , down-weighting each over-sampled condition"%(str(len(all_evoked_condition_A)),str(len(all_evoked_condition_B))))
            all_evoked_condition_A = all_evoked_condition_A[0:min_dim]
            all_evoked_condition_B = all_evoked_condition_B[0:min_dim]            


        diff_evokeds = [mne.combine_evoked([all_evoked_condition_B[r],all_evoked_condition_A[r]],weights=[1,-1]) for r in range(0,len(all_evoked_condition_A))] 
        diff_evokeds_group = mne.combine_evoked([diff_evoked for diff_evoked in diff_evokeds], weights='equal') #[1 for diff_evoked in diff_evokeds] );
        evo_A = mne.combine_evoked(all_evoked_condition_A,weights='equal')
        evo_B = mne.combine_evoked(all_evoked_condition_B,weights='equal')
        diff_evokeds_group = mne.combine_evoked([evo_B,evo_A], weights=[1,-1])

        

        #### MAIN CHOICE AS TO WHAT TO PLOT
        #headshape = diff_evokeds_group
        headshape =  copy.deepcopy(self.group_A_avg['B'])    # Young person's head      
        #headshape = centre_sensor_locations(headshape)[0]
        info = headshape.info
        pos = mne.find_layout(info).pos
        pos_l = list(pos)
        # Delete location
        chan_nums_to_delete = sorted([int(x.replace("MEG ","")) for x in self.missing_channels],reverse=True)
        print("Chans to delete: ", chan_nums_to_delete)        
        #pos_l[(chan_nums_to_delete[0]-1):(chan_nums_to_delete[0]+1)]
        for index in chan_nums_to_delete:
            del pos_l[index]
        pos_l = pos_l[0:NUM_CHANS_TO_KEEP]        
        # After
        pos = np.array(pos_l)
        # Now it's theoretically including bad channels....
        print("Number of channels not bad :", len(pos))

        # configure variables for visualization
        colors = {"rel_target": "crimson", "unrel_target": 'steelblue'}
        # get sensor positions via layout
        # find offset centroid and remove
        x = [p[0] for p in pos[:,0:2]]
        y = [p[1] for p in pos[:,0:2]]
        centroid = (sum(x) / len(pos), sum(y) / len(pos))
        # Scale by a factor of 5 (this is weird)
        pos_scaled_centred = np.array([[(points[0]-centroid[0])/5,(points[1]-centroid[1])/5,points[2],points[3]] for points in pos])
        
        # Adjacency matrix
        ch_adjacency, ch_names = mne.channels.find_ch_adjacency(info, 'mag');

        

        ## Structure the data to do the test over
        X = np.array([k.data for k in all_evoked_condition_A]), np.array([j.data for j in all_evoked_condition_B])  # as 3D matrix
        X = [np.transpose(x, (0, 2, 1)) for x in X]  # transpose for clustering
        

        # Set family-wise p-value
        t_obs, clusters, cluster_pv, h0 = mne.stats.spatio_temporal_cluster_test(X, n_permutations=N_PERMUTATIONS, 
                                                     #threshold=threshold, # Inputting NONE makes it use p < 0.05
                                                     tail=0,
                                                     n_jobs=1,
                                                     adjacency=ch_adjacency
                                                     ) 
        print("Cluster p values ", cluster_pv)
        good_cluster_inds = np.where(cluster_pv < CLUSTER_CUTOFF)[0]
        print("Good cluster indices: ", good_cluster_inds)






        if mode == 'condition':
            extra_filename = 'Cond B - Cond A '
        if mode == "mmr_by_age":
            extra_filename = "Ages max "+str(age_cutoff_low) + " min "+str(age_cutoff_high)+" ("+str(size_age_bucket_percent)+"%)"
            extra_filename +=" MMR vs age "              
        elif "age" in mode:
            extra_filename = "Ages max "+str(age_cutoff_low) + " min "+str(age_cutoff_high)+" ("+str(size_age_bucket_percent)+"%)"
            extra_filename += " Condition "
            if condition_comparison == None:
                extra_filename += "EveryTrial "      
            else:
                extra_filename += str(condition_comparison)+" "  
        xlabel = ""
        if mode == "condition":
            xlabel+= "Mean difference Condition B minus Condition A\n"
        if mode in ["age", "condition_by_age","mmr_by_age"]:
            xlabel+= "Mean difference Oldest minus youngest "+str(size_age_bucket_percent)+"%,\n"
        if mode in ["age", "condition_by_age"]:            
            xlabel+=str("Condition: ")
            if condition_comparison!=None:
                xlabel+=str(condition_comparison)+" "
            else:
                xlabel+=str("EVERY TRIAL ")
        if mode in ["mmr_by_age"]:
            xlabel+=str("in the MMR: ")

        filename_all_clusters = 'SpatioTemporalCluster ERF #All '+str(extra_filename)+str(int(self.participants[0].epochs_ransac_autoreject.times[0]*SEC_TO_MS+EPSILON_TIME))+ "-"+str(int(self.participants[0].epochs_ransac_autoreject.times[-1]*SEC_TO_MS+EPSILON_TIME))+"ms Details"
        with open(filename_all_clusters+'.txt', 'w') as g:
            g.write("Good cluster indices\n")
            g.write(str(good_cluster_inds))
            g.write("\n")
            g.write("All cluster p values\n")
            g.write(str(cluster_pv))



        times = headshape.to_data_frame()['time'].values #  diff_evokeds_group.to_data_frame()['time'].values
        for i_clu, clu_idx in enumerate(good_cluster_inds):
            time_inds, space_inds = np.squeeze(clusters[clu_idx])
            ch_inds = np.unique(space_inds)
            time_inds = np.unique(time_inds) 
            # find the time points of significance
            print("Significant times", [times[r] for r in time_inds]) #[sample_time+5*time_ind for time_ind in time_inds])        
            sig_times = self.participants[0].epochs_ransac_autoreject.times[time_inds]
            #print("Significant times", sig_times)
            sig_chans = ["MEG "+str(ch_ind+1).zfill(3) for ch_ind in ch_inds]                
            print("Significant channels ", sig_chans )

            # get topography for F stat
            t_map = t_obs[time_inds, ...].mean(axis=0)

            # get topography of difference
            time_shift = self.participants[0].epochs_ransac_autoreject.time_as_index(self.participants[0].epochs_ransac_autoreject.times[0])      # fix windowing shift
            
            #diff_topo = np.mean(headshape.data[:,time_inds+time_shift],axis=1)
            #diff_topo = np.mean(diff_evokeds.data[:,time_inds+time_shift],axis=1)
            diff_topo = np.mean(diff_evokeds_group.data[:,time_inds+time_shift],axis=1)        
            

            # create spatial mask
            mask = np.zeros((t_map.shape[0], 1), dtype=bool)
            mask[ch_inds, :] = True

            # initialize figure
            fig, ax_topo = plt.subplots(nrows=1, ncols=1, figsize=(20, 6)) # ncols=2,
            fig,(ax1,ax2) = plt.subplots(ncols=2)

            #shifts = centre_sensor_locations(headshape) # centre_sensor_locations(diff_evokeds_group)
            #shift_x,shift_y,shift_z = shifts[1], shifts[2], 0

            # Plot average difference and mark significant sensors
            image, _ = plot_topomap(diff_topo,
                                    pos_scaled_centred,
                                    mask=mask, axes=ax_topo, cmap='RdBu_r',
                                    #title=None,
                                    #sphere=(shift_x, shift_y, shift_z, CONSTANT_RADIUS_MULT*find_radius(diff_evokeds_group)), # find_radius(diff_evokeds_group)), 
                                    res=RESOLUTION,
                                    vmin=-1*10**(-14), vmax=1*10**(-14),
                                    show=False,
                                    #vmin=np.min, vmax=np.max, show=False,
                                    outlines='head')

            # Create additional axes (for ERF and colorbar)
            divider = make_axes_locatable(ax_topo)

            # Add axes for colorbar
            ax_colorbar = divider.append_axes('right', size='10%', pad=0.05)
            plt.colorbar(image, cax=ax_colorbar)

            ax_topo.set_xlabel(
               xlabel+'({:0.3f} - {:0.3f} s)'.format(*sig_times[[0, -1]]))


        # Add new axis for time courses and plot time courses
            ax_signals = divider.append_axes('right', size='300%', pad=1.2)
            #ax_signals.set_xlim([PLOT_START[0]*1/SEC_TO_MS, times[-1]*1/SEC_TO_MS+EPSILON_TIME])
            #ax_signals.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS, xmax=EPOCH_END_ANALYSIS*SEC_TO_MS)
            ax_signals.set_xlim(xmin=PLOT_START*SEC_TO_MS, xmax=PLOT_END*SEC_TO_MS)
            title = 'Cluster #{0}, {1} sensor'.format(i_clu + 1, len(ch_inds))
            if len(ch_inds) > 1:
                title+="s"
            title+=" had significant differences"


            plot_compare_evokeds(diff_evokeds_group, # evokeds
                            title=title,
                            picks=ch_inds,
                            axes=ax_signals, # 
                            #colors=colors,
                            show=False,
                            combine='mean',
                            legend=True,
                            #show_sensors=True,
                            ci=True,
                            #split_legend=True,
                            #truncate_yaxis='max_ticks'
                           )



            # plot temporal cluster extent
            ymin, ymax = ax_signals.get_ylim()
            ax_signals.fill_betweenx((ymin, ymax), sig_times[0], sig_times[-1],
                                    color='orange', alpha=0.3)

            # clean up viz
            mne.viz.tight_layout(fig=fig)
            fig.subplots_adjust(bottom=.05)
            plt.legend('',frameon=False)      

            try:
                ax_signals.get_legend().remove()
            except Exception as e:
                print("Cannot remove legend ", e)


            # save the plot as a file
            filename = 'SpatioTemporalCluster ERF #'+str(clu_idx+1)+' '+str(extra_filename)+str(int(self.participants[0].epochs_ransac_autoreject.times[0]*SEC_TO_MS+EPSILON_TIME))+ "-"+str(int(self.participants[0].epochs_ransac_autoreject.times[-1]*SEC_TO_MS+EPSILON_TIME))+"ms"
            #ax_signals.figure
            ax_signals.figure.savefig(filename+'.jpg',
                        format='jpeg',
                        dpi=DPI,
                        bbox_inches='tight')

            with open(filename+' Details.txt', 'w') as f:
                f.write("Chans\n")
                f.write(str(sig_chans)+"\n")
                f.write("Times\n")
                f.write(str([times[r] for r in time_inds]))
                f.write("\n")
                f.write("\n")                
                f.write("t-map\n")
                f.write(str(t_map))
                f.write("\n")
                f.write("t-map at sig channels\n")
                try:
                    for ch in sig_chans:
                        #print("Ch ", ch)
                        ch_ = ch.replace("MEG ","")
                        ch_ = int(ch_)
                        ch_index = ch_ - 1
                        f.write(str(t_map[ch_index])+", ")
                except:
                    self.warn("ERROR WRITING DETAILS"+str(sig_chans))
                f.write("\n")
                f.write("\n")           
                f.write("t-obs\n")
                f.write(str(t_obs))
                f.write("\n")
                f.write("\n")               
                f.write("Mask\n")                
                f.write(str(mask))





    def statistical_head_spatio_temporal(self,method_dct,waveform,age_group="Both",stat_to_use=STAT_TO_USE):
        
        '''
        Generates topo plots of spatiotemporal clusters where a statistic (such as correlation) of the regression between predictor and ERF is significant
        Also generates evoked graphs during these times of significant difference
        Also outputs text files stating on which channels the difference is significant
        '''
        
        # Doing stats on all age groups?
        if age_group=="Both":
            file_name_addition=""
        else:
            file_name_addition=" "
            if age_group=="Young":
                file_name_addition+= '%.2f'%(self.age_young_min)+'-'+'%.2f'%(self.age_young_max)+" years"
            elif age_group=="Old":
                file_name_addition+= '%.2f'%(self.age_old_min)+'-'+'%.2f'%(self.age_old_max)+" years"
            





        print("@@@@@@@@@@@@@@ OBTAINING STATISTICS TOPO PLOTS @@@@@@@@@@@@@@")
        self.remove_surplus_channels()
        self.identify_missing_channels()
        self.get_child_adult_sample()
        condition_A = self.sample_child.cond_A
        condition_B = self.sample_child.cond_B

        # Run analysis getting correlations only so we can test these for significant differences to zero
        # Prepare data structures
        stats            = {}
        gfps_all_ptcps   = {}
        stats_all_ptcps  = {}
        ptcps_surprise   = {}
        all_predictions  =  []
        all_predictions_chan  = {}
        all_predictions_chan_all_ptcps = {}
        all_predictions_total = {}
        all_values_chan       = {}
        all_values_gfps       = {}
        value_dct_all_ptcps   = {}
        for statistic in [stat_to_use]:
            stats_all_ptcps[statistic] = {}
            stats[statistic] = {}   
        
        exception_log    = []
        length_divisions = []  # Keep track of how many epochs we use for each participant in this analysis
        j = 0
        k = 0
        cum_channels = 0
        total_number_events = 0
        sample_time = int(EPOCH_START_ANALYSIS*SEC_TO_MS)
        for ptcp in self.participants:

            include_in_analysis = False
            if age_group=="Both":
                include_in_analysis=True
            elif age_group=="Young":
                if ptcp.p_id in self.youngest:
                    include_in_analysis=True
            elif age_group=="Old":
                if ptcp.p_id in self.oldest:
                    include_in_analysis=True
            else:
                self.warn("YOU HAVE NOT SPECIFIED AGE CUTOFF FOR statistical_head_spatio_temporal ")

            if include_in_analysis:
                try:
                    print("########### ", ptcp.p_id, "Ptcp Included #: ", j)


                    epochs_all = ptcp.epochs_ransac_autoreject
                    epochs_all_df = epochs_all.to_data_frame()
                    L = len(epochs_all)    

                    ## Get times
                    times = [int(item*SEC_TO_MS+EPSILON_TIME) for item in ptcp.epochs_ransac_autoreject.times]
                    #print("Times are ", len(times))

                    ptcp.get_surprise_hgf()    # Unnecessary to reload these, should already have it
                    ptcp.surprise_prediction = ptcp.surprise[ptcp.surprise_relevant_series]
                    pred_surprise_vals = ptcp.surprise_prediction    # Uses self.surprise_prediction = self.surprise[self.surprise_relevant_series]

                    if method_dct['On every trial'] == True:
                        epochs_all_df = epochs_all.to_data_frame()
                        epochs_all_surprise_indices = sorted(list(set(epochs_all_df['epoch'].values)))
                        retained_epoch_indices = epochs_all_surprise_indices # retain all trials
                    elif method_dct['On low surprise trials'] == True:
                        epochs_low_surprise = epochs_all[condition_A]
                        epochs_low_surprise_df = epochs_low_surprise.to_data_frame()
                        epochs_low_surprise_indices = sorted(list(set(epochs_low_surprise_df['epoch'].values)))
                        retained_epoch_indices = epochs_low_surprise_indices
                        pred_surprise_vals = [pred_surprise_vals[r] for r in epochs_low_surprise_indices] # Subset of predictions          
                    elif method_dct['On high surprise trials'] == True:       
                        epochs_high_surprise = epochs_all[condition_B]
                        epochs_high_surprise_df = epochs_high_surprise.to_data_frame()
                        epochs_high_surprise_indices = sorted(list(set(epochs_high_surprise_df['epoch'].values)))
                        retained_epoch_indices = epochs_high_surprise_indices
                        pred_surprise_vals = [pred_surprise_vals[r] for r in epochs_high_surprise_indices] # Subset of predictions              
                    elif method_dct['Subtracting deviants only'] == True or method_dct['Subtracting high surprise average'] == True :     
                        # Method A to calculate which epochs to retain
                        evts = ptcp.events[:, 2]
                        evts_high_surprise = [r for r in range(0,len(evts)) if evts[r] == float(condition_B)]
                        # pred_surprise_vals = [pred_surprise_vals[r] for r in evts_high_surprise] # Subset of predictions
                        # retained_epoch_indices = evts_high_surprise

                        # Method B to calculate which epochs to retain
                        epochs_high_surprise = epochs_all[condition_B]
                        epochs_high_surprise_df = epochs_high_surprise.to_data_frame()
                        epochs_high_surprise_indices = sorted(list(set(epochs_high_surprise_df['epoch'].values)))
                        retained_epoch_indices = epochs_high_surprise_indices
                        pred_surprise_vals = [pred_surprise_vals[r] for r in epochs_high_surprise_indices] # Subset of predictions
                    elif method_dct['Subtracting deviants only'] == True:                   
                        # Method A to calculate which epochs to retain
                        evts = ptcp.events[:, 2]
                        evts_low_or_high = [r for r in range(0,len(evts)) if (evts[r] == float(condition_A) or evts[r] == float(condition_B))]


                        # Method B to calculate which epochs to retain
                        epochs_low_or_high = epochs_all[condition_B]
                        epochs_low_or_high_df = epochs_low_or_high.to_data_frame()
                        epochs_low_or_high_indices = sorted(list(set(epochs_low_or_high_df['epoch'].values)))
                        retained_epoch_indices = epochs_low_or_high_indices
                        pred_surprise_vals = [pred_surprise_vals[r] for r in epochs_low_or_high_indices] # Subset of predictions
                    elif method_dct['Subtracting low surprise average'] == True: 
                        # Method A to calculate which epochs to retain
                        evts = ptcp.events[:, 2]
                        evts_low_surprise = [r for r in range(0,len(evts)) if evts[r] == float(condition_A)]
                        # pred_surprise_vals = [pred_surprise_vals[r] for r in evts_low_surprise] # Subset of predictions
                        # retained_epoch_indices = evts_low_surprise

                        # Method B to calculate which epochs to retain
                        epochs_low_surprise = epochs_all[condition_A]
                        epochs_low_surprise_df = epochs_low_surprise.to_data_frame()
                        epochs_low_surprise_indices = sorted(list(set(epochs_low_surprise_df['epoch'].values)))
                        retained_epoch_indices = epochs_low_surprise_indices
                        pred_surprise_vals = [pred_surprise_vals[r] for r in epochs_low_surprise_indices] # Subset of predictions
                    else:
                        retained_epoch_indices = list(set(epochs_all_df['epoch'].values))
                        pred_surprise_vals = [pred_surprise_vals[i] for i in retained_epoch_indices] # The predictions pertaining to the epochs retained
                    L = len(retained_epoch_indices)
                    print("Number of events analysing: ", L )
                    length_divisions.append(L)
                    total_number_events = total_number_events  + L


                    #Different otpions for which signal to correlate with the predicted surprise value
                    ##Option 1 - all epochs, unless over-written
                    ptcp.epochs_ransac_autoreject.baseline = (EPOCH_START_ANALYSIS,EPOCH_END_ANALYSIS)
                    conditionA = ptcp.epochs_ransac_autoreject[condition_A]
                    conditionB = ptcp.epochs_ransac_autoreject[condition_B]
                    # print("Con A", len(conditionA.to_data_frame().index.values))
                    # print("Con B", len(conditionB.to_data_frame().index.values))            
                    concat_epochs = mne.concatenate_epochs([conditionA,conditionB])          
                    conditionAll = concat_epochs
                    # print("Con All", len(conditionAll.to_data_frame().index.values))         
                    
                    # print(len(conditionA),len(conditionB))
                    ## Option 2 - evoked with code 
                    if method_dct['On every trial'] == True:
                        epochs_diff = conditionAll
                    elif method_dct['On low surprise trials'] == True:
                        epochs_diff = conditionA            
                    elif method_dct['On high surprise trials'] == True:
                        epochs_diff = conditionB            
                    elif method_dct['Subtracting average of all trials'] == True:
                        # The epoch object
                        epochs_df  = ptcp.evoked_generic_dfs  # All
                        #epochs_df = ptcp.evoked_generic_dfs[cond_code] # doesn't work as this is an average
                        all_epochs = ptcp.epochs_ransac_autoreject
                        all_epochs_evoked = all_epochs.average() 

                        epochs_diff = all_epochs.subtract_evoked(all_epochs_evoked)     # Subtract avg from each epoch
                    elif method_dct['Subtracting low surprise average'] == True:
                        conditionA_evoked = conditionA.average()                

                        epochs_diff = conditionA.subtract_evoked(conditionA_evoked)     # Subtract the average of all condition A (pre-deviants) from condition A(expected to be unsurprising) for each epoch                
                    elif method_dct['Subtracting high surprise average'] == True:
                        conditionB_evoked = conditionB.average()                

                        epochs_diff = conditionB.subtract_evoked(conditionB_evoked)     # Subtract the average of all condition B (deviants) from condition B(expected surprises) for each epoch
                    elif method_dct['Subtracting deviants only'] == True:
                        conditionA_evoked = conditionA.average() 
                        conditionB_evoked = conditionB.average()

                        epochs_diff = conditionB.subtract_evoked(conditionA_evoked)     # Subtract condition A(expected to be unsurprising) from condition B(expected surprises) for each epoch


                    epochs_df   = epochs_diff.to_data_frame()        
                    # print("Epochs df length ", len(epochs_df))
                    # Replace strings as necessary
                    try:
                        if type(pred_surprise_vals[0])==type("0.5"):
                            pred_surprise_vals_exp = [float(x.replace("\n","")) for x in pred_surprise_vals]
                        elif type(pred_surprise_vals[0])==type(0.5):
                            pred_surprise_vals_exp = pred_surprise_vals
                    except Exception as e:
                        print("PROBLEM GENERATING PREDICTED SURPRISE VALUES", str(e))      
                    ptcps_surprise[ptcp.p_id] = pred_surprise_vals_exp
                    #print("Len surprise ", len(pred_surprise_vals_exp))
                    value_dct   = {}
                    gfps        = {}
                    for statistic in [stat_to_use]: # Prepare data structures
                        stats[statistic] = {}    
                        stats_all_ptcps[statistic][ptcp.p_id] = {}


                    for time_ in times:
                        #print("Time ", time_)
                        if time_ not in list(all_predictions_chan.keys()): # Prepare data structures
                            print(time_, " Adding time to dct")
                            all_predictions_chan[time_]  = []
                            all_predictions_total[time_] = []
                            all_values_chan[time_]       = []
                            all_values_gfps[time_]       = []
                        if time_ not in list(all_predictions_chan_all_ptcps.keys()):
                            all_predictions_chan_all_ptcps[time_] = []

                        if time_ >= EPOCH_START_ANALYSIS*SEC_TO_MS and time_ <= EPOCH_END_ANALYSIS*SEC_TO_MS:
                            if time_ % 30 == 0: # Progress update
                                print("Time : ", time_)
                            for statistic in [stat_to_use]:
                                stats[statistic][time_] = {}
                                stats_all_ptcps[statistic][ptcp.p_id][time_] = {}

                            epochs_time_snapshot = epochs_df[epochs_df['time'] == time_]  # Multiple snapshots at the same time in different epochs
                            #print("Length of epochs time snapshot at time %s is %s "%(str(time_),str(len(epochs_time_snapshot))))
                            
                            #Every column individually
                            ttl = 0
                            num_chans = 0
                            for chan in epochs_time_snapshot.columns.values:
                                if "MEG" in chan:
                                    chan_values = epochs_time_snapshot[chan].values   # list of 
                                    ttl+=len(chan_values)
                                    num_chans+=1
                                    if chan in value_dct.keys():
                                        if time_ in value_dct[chan]:
                                            value_dct[chan][time_].append(chan_values)
                                        else:
                                            value_dct[chan][time_] = chan_values
                                    else:
                                        value_dct[chan] = {}
                                        value_dct[chan][time_] = chan_values
                                    if chan in value_dct_all_ptcps.keys(): 
                                        if time_ in value_dct_all_ptcps[chan].keys():
                                            value_dct_all_ptcps[chan][time_] = list(value_dct_all_ptcps[chan][time_])
                                            for chan_val in chan_values:
                                                value_dct_all_ptcps[chan][time_].append(chan_val)
                                            value_dct_all_ptcps[chan][time_] = np.array(value_dct_all_ptcps[chan][time_])
                                            #print("In, length afterwards ", len(value_dct_all_ptcps[chan][time_]))
                                        else:
                                            value_dct_all_ptcps[chan][time_] = chan_values
                                            #print("Not in, length afterwards ", len(value_dct_all_ptcps[chan][time_]))
                                    else:
                                        value_dct_all_ptcps[chan] = {}
                                        value_dct_all_ptcps[chan][time_] = chan_values

                                    all_predictions_chan_all_ptcps[time_].append(pred_surprise_vals_exp)     # The prediction is the same for every channel
                                    all_predictions_chan[time_].append(pred_surprise_vals_exp)
                                    all_values_chan[time_].append(value_dct[chan][time_])


                                    # Available statistics
                                    corr = pearsonr(pred_surprise_vals_exp, value_dct[chan][time_])[0]
                                    r_squared = corr**2
                                    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(pred_surprise_vals_exp, value_dct[chan][time_])
                                    if len(pred_surprise_vals_exp) != len(chan_values):
                                        self.warn("Length mismatch " +str(time_) + str(" ")+str(chan) +str(" ")+ str(len(pred_surprise_vals_exp)) + str(" ") + str(len(chan_values) ))

                                    # Save all computed statistics
                                    for statistic in [stat_to_use]:
                                        stat = eval(statistic) # corr # r_squared # p_value # slope            
                                        stats[statistic][time_][chan] = stat # Only works for one ptcp
                                        stats_all_ptcps[statistic][ptcp.p_id][time_][chan] = stat


                            # Sum over columns    
                            cols_to_keep = []
                            for col in epochs_time_snapshot.columns.values:
                                if "MEG" in col:
                                    cols_to_keep.append(col)
                            epochs_time_snapshot = epochs_time_snapshot[cols_to_keep]

                            gfps[time_] = []       
                            for idx in range(0,len(epochs_time_snapshot)):#.index.values:
                                chan_vals = epochs_time_snapshot.iloc[idx]
                                square_array = [chan**2 for chan in chan_vals.values]
                                total = sum(square_array)
                                mean_sq = total/len(square_array)
                                mean_sq_sqrt = math.sqrt(mean_sq)
                                gfps[time_].append(mean_sq_sqrt)
                                all_predictions_total[time_].append(pred_surprise_vals_exp[idx])
                                all_values_gfps[time_].append(mean_sq_sqrt)        

                    
                    gfps_all_ptcps[ptcp.p_id] = gfps
                    ptcp_predictions = all_predictions_chan_all_ptcps[sample_time][cum_channels]   # # The prediction is the same for every channel and time (in the epoch), so e.g. all_predictions_chan_all_ptcps[sample_time][0] and all_predictions_chan_all_ptcps[sample_time][124]  are identical    
                    for pred in ptcp_predictions:
                        all_predictions.append(pred)  # Pick an arbitrary time (here, 130ms) and channel (here, "MEG 005") as the prediction is the same at every channel
                    print("Length of all predictions ", len(all_predictions))
                    assert(len(all_predictions)==len(value_dct_all_ptcps['MEG 011'][sample_time])), "Mismatch in lengths"
                    k+=1
                    cum_channels+=num_chans
                    #clear_output()
                    print("Exception log ", exception_log)


                except Exception as e:
                    print("@@@@@@ PROBLEM WITH PARTICIPANT ", ptcp.p_id, e)
                    exception_log.append([ptcp.p_id, e])
                    #clear_output()
                    for exception_logged in exception_log:
                        print(exception_logged)
                j+=1

        try_del(['epochs_all_showing_times', 'epochs_all_showing_times_df', 'epochs_time_snapshot', 'epochs_all_df', 'total', 'totals'])
        gc.collect()




        ### NOTE THIS IS FOR CONDITION B only
        for stat in [stat_to_use]: 
            # Intermediate step
            stat_arrays_int = {}
            count = 0
            for ptcp in stats_all_ptcps[stat].keys():
                print("Times, ", stats_all_ptcps[stat][ptcp].keys())
                for time in stats_all_ptcps[stat][ptcp].keys():
                    for chan in stats_all_ptcps[stat][ptcp][time].keys():
                        if chan not in stat_arrays_int.keys():
                            stat_arrays_int[chan] = {}  

                        stat_val = stats_all_ptcps[stat][ptcp][time][chan]
                        stat_val = np.arctanh(stat_val) if stat == "corr" else stat_val
                        if time not in stat_arrays_int[chan].keys():
                            stat_arrays_int[chan][time] = [stat_val]  
                        else:
                            stat_arrays_int[chan][time].append(stat_val)

            for chan in stat_arrays_int.keys():
                for time in stat_arrays_int[chan].keys():
                    avg_stat = sum(stat_arrays_int[chan][time])/len(stat_arrays_int[chan][time])
                    if stat == corr:
                        avg_stat = np.tanh(avg_stat)  # /(k+1)      
                    stat_arrays_int[chan][time] = avg_stat

            # Lay out these in the format required to plot topographically
            stats_arrays = {}
            #print("Stat arrays int ", stat_arrays_int)
            first_key = list(stat_arrays_int.keys())[0]
            times = stat_arrays_int[first_key].keys()
            for chan in stat_arrays_int.keys():
                stats_arrays[chan] = [stat_arrays_int[chan][time] for time in times]


            # Put this onto a headshape
            array_avg_stat = []
            for chan in list(stats_arrays.keys()):
                array_avg_stat.append([stats_arrays[chan][r] for r in range(0,len(sorted(times))) ])
            self.get_child_adult_sample()
            statistics_head_avg_stat = self.sample_child.epochs_ransac_autoreject[self.sample_child.cond_B].average() # Young people smaller head
            #statistics_head_avg_stat = self.sample_child.epochs_ransac_autoreject[self.sample_child.cond_B]  # Young people smaller head
            #statistics_head_avg_stat = copy.deepcopy(self.group_A_avg['B'])                             # Young people smaller head
            array_num = 0
            for array in array_avg_stat:
                array = np.array([x* 10**(-15) for x in array])         # Was in fT
                array_avg_stat[array_num] = array
                array_num+=1
            statistics_head_avg_stat.data = np.array(array_avg_stat ) 


            title = ' averaged ' +str(stat) + ' over all participants between predictor and ERF, ' + waveform
            cmap = 'RdBu_r' if stat!='p_value' else 'RdBu' # Reverse colour scheme for p_values (lower is more significant, which is represented by red colour)      # -> NEED TO CORRECT FOR MULTIPLE COMPARISONS

            #shifts = centre_sensor_locations(statistics_head_avg_stat)
            #shift_x,shift_y,shift_z = shifts[1], shifts[2], 0

        for statistic in [stat_to_use]:
            data = stats_all_ptcps[statistic]
            for ptcp in data.keys():
                data_ptcp = data[ptcp]
                

        REF_OBJECT = mne.combine_evoked([self.group_B_avg['B'],self.group_B_avg['A'], self.group_A_avg['B'],self.group_A_avg['A']], weights=[-1,1,1,-1])
        chans_data = REF_OBJECT.data
        times_data = REF_OBJECT.data[0]
                
        N = len(chans_data)
        M = len(times_data)
        stat_on_head = [[] for x in range(0,N)]         # [np.array([]) for x in range(0,N)]   
        #### NOTE THIS IS FOR CONDITION B ONLY
        for stat in [stat_to_use]:
            stat_on_head = [[] for x in range(0,N)]     # [np.array([]) for x in range(0,N)]

            # Construct head object
            m = 0
            print("Stat :", stat)
            for chan in self.group_A_avg['B'].info['ch_names']: # 
                #print("Chan ", chan)
                stat_on_head[m] = []
                for time in [int(item*SEC_TO_MS+EPSILON_TIME) for item in self.participants[-1].epochs_ransac_autoreject.times]: #stats_all_ptcps[stat][self.sample_adult.p_id].keys(): # times = [int(item*SEC_TO_MS+EPSILON_TIME) for item in ptcp.epochs_ransac_autoreject.times]

                    corr = pearsonr(all_predictions, value_dct_all_ptcps[chan][time])[0]  # value_dct_all_ptcps)[chan][time]      
                    r_squared = corr**2
                    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(all_predictions,  value_dct_all_ptcps[chan][time]) # value_dct_all_ptcps[chan][time])
                    relevant_stat = eval(stat)
                    stat_on_head[m].append(relevant_stat)
                m+=1


            # Replace the channel ERF data at any time/chan with the correlation for that time/chan instead of ERF
            self.get_child_adult_sample()
            statistics_head     = copy.deepcopy(self.sample_child.epochs_ransac_autoreject[self.sample_child.cond_A].average()) # Young people smaller head
            #statistics_head    = copy.deepcopy(self.group_A_avg['B'])                                                  # Young people smaller head
            # statistics_head   = centre_sensor_locations(statistics_head)[0]

            array_num = 0
            for array in stat_on_head:
                #array = np.array([abs(x) for x in array])
                array = np.array([x* 10**(-15) for x in array])        
                stat_on_head[array_num] = array
                array_num+=1
            statistics_head.data = np.array(stat_on_head) #  * 10**(-15) # Was in fT

            title = stat + ' between predictor and ERF, ' + waveform + ', for whole group data'
            if stat == "p_value":
                title = 'P-value of linear regression slope between predictor and ERF at various times, for whole group data'
            if stat == "slope":
                title = 'Slope of linear regression between predictor and ERF at various times, for whole group data'

            cmap = 'RdBu_r' if stat!='p_value' else 'RdBu' # Reverse colour scheme for p_values (lower is more significant, which is represented by red colour)      # -> NEED TO CORRECT FOR MULTIPLE COMPARISONS
            #cmap = 'RdBu_r'

            vmin = 0 if stat == 'p_value' else None
            vmax = 1 if stat == 'p_value' else None


            #shifts = centre_sensor_locations(statistics_head)
            #shift_x,shift_y,shift_z = shifts[1], shifts[2], 0


        ## SINGLE STATISTIC
        # Correlation across all participants at each time point, aggreating over all channels for each epoch
        first_key = list(gfps_all_ptcps.keys())[0]
        all_times = list(gfps_all_ptcps[first_key].keys())
        preds = {}
        signals = {}
        corr_saved_values = []
        slope_saved_values = []
        t_saved_values = []
        z_saved_values = []
        z_confidences = []
        r_confidences = []
        for time in all_times:
            preds[time] = []
            signals[time] = []
            for j in range(0,len(all_values_gfps[time])):
                predictions = [all_values_gfps[time][j]]        # Gfps at each time
                meg_signal = [all_predictions_total[time][j]]   # Predicted surprise value at each time (the same for every time)
                for pred in predictions:
                    preds[time].append(pred)
                for signal in meg_signal:
                    signals[time].append(signal)     
            corr = pearsonr(preds[time], signals[time])[0]
            corr_saved_values.append(corr)
            # Other stats
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(preds[time], signals[time])
            slope_saved_values.append(slope)

            n = len(preds[time])
            print("Number of samples for correlation confidence interval is ", n)
            t_stat = corr_coeff_t_value(n=n,corr=corr)        
            t_saved_values.append(t_stat)

            corr_null = 0
            #z_observed = 0.5*math.log((1+corr)/(1-corr)) # fisher r to z transformation
            z_observed = fisher_r_to_z(corr)
            z_null = 0.5*math.log((1+corr_null)/(1-corr_null))
            z = (z_observed-z_null) / math.sqrt(1/(n-3))
            z_saved_values.append(z)     
            
            bonferroni_adjusted_cutoff = CLUSTER_CUTOFF/2 / (4*13)
            z_critical = st.norm.ppf(1 - bonferroni_adjusted_cutoff) # 1.96 at p=0.05
            z_confidence = [z_null-z_critical*math.sqrt(1/(n-3)),z_null+z_critical*math.sqrt(1/(n-3))]
            r_confidence = [fisher_z_to_r(z_confidence[0]), fisher_z_to_r(z_confidence[1])]

            print("Time: ",  time, " Correlation: ", corr, "z_obs", z_observed, " z score ", z)
            # print("t saved values ", t_saved_values)
            # print("z saved values ", z_saved_values)
            # print("Z confidence ", z_confidence)
            print("r confidence ", r_confidence)
            z_confidences.append(z_confidence)
            r_confidences.append(r_confidence)

        from matplotlib import pyplot as plt
        try:
            plt.close('all')
            plt.clf()
        except:
            print("Problem closing graphs")
        if stat_to_use == 'corr':
            return_val = corr_saved_values
            title = "Correlation of predictor to GFP\n"+waveform+file_name_addition
        if stat_to_use == 'slope':
            return_val = slope_saved_values 
            title = "Slope of regression (predictor to GFP)\n"+waveform+file_name_addition
        ylab = "Correlation" if stat_to_use == 'corr' else "Slope"


        plt = plot_module(  xs = [all_times],
                            ys = [return_val],
                            xlabel = "Time (ms)", ylabel=ylab,
                            legends=[stat_to_use],
                            title=title,
                            font_object = self.font,
                            x_min_to_show = PLOT_START*SEC_TO_MS,
                            x_max_to_show = PLOT_END*SEC_TO_MS,                              
                            #x_min_to_show = EPOCH_START_ANALYSIS*SEC_TO_MS,
                            #x_max_to_show = EPOCH_END_ANALYSIS*SEC_TO_MS                            
                            )
        # plt.plot(all_times,return_val)
        # plt.title(title, fontsize=FONTSIZE_TITLE)        
        # plt.ylabel(ylab, fontsize=FONTSIZE_LABELS) 
        # plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
        # #plt.xlim([min(all_times),max(all_times)])
        # plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed 

        timestring = str(int(self.participants[0].epochs_ransac_autoreject.times[0]*SEC_TO_MS+EPSILON_TIME))+ "-"+str(int(self.participants[0].epochs_ransac_autoreject.times[-1]*SEC_TO_MS+EPSILON_TIME))+"ms"
        filename = 'Relationship of '+ylab +' to GFP, '+ waveform + ' ' + timestring+' '+file_name_addition
        plt.savefig(filename+'.jpg',
                    format='jpeg',
                    dpi=DPI,
                    bbox_inches='tight') 
        plt.show()  




    #     plt.plot(all_times,t_saved_values)
    #     plt.title("t value of corr HGF_Vanilla PE2 to GFP, 20pct threshold, "+waveform, fontsize=FONTSIZE_TITLE)
    #     plt.ylabel("t statistic", fontsize=FONTSIZE_LABELS) 
    #     plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
    #     plt.show()    
        


        if age_group == 'Old':
            sample_individual = self.sample_adult
        else:
            sample_individual = self.sample_child

        # Build an array with the right dimensions
        num_participants = len(stats_all_ptcps[stat_to_use].keys())
        ptcp_array = [np.array([]) for x in range(0,num_participants)]  
        num_times = len(stats_all_ptcps[stat_to_use][sample_individual.p_id].keys())
        times_array = [np.array([]) for x in range(0,num_times)]  
        num_chans = len(stats_all_ptcps[stat_to_use][sample_individual.p_id][sample_time])
        chans_array = [np.array([]) for x in range(0,num_chans)]
        new_X = [np.array([]) for x in range(0,num_participants)]
        p = 0
        for ptcp_array in new_X:
            new_X[p] = times_array
            p+=1

        p = 0
        for ptcp_array in new_X:
            c = 0
            for time in new_X[p]:
                new_X[p][c] = chans_array
                c+=1
            p+=1    
        print("New X shape", np.shape(new_X))

        # Populate the values
        p = 0
        for ptcp in stats_all_ptcps[stat_to_use].keys(): # These do not need to be arranged by time/space, it's just participant by participant
            t = 0
            for time in sorted(stats_all_ptcps[stat_to_use][ptcp]): # These are temporally arranged
                c = 0
                for chan in stats_all_ptcps[stat_to_use][ptcp][time].keys(): # These are not spatially arranged
                    new_X[p][t][c] = stats_all_ptcps[stat_to_use][ptcp][time][chan]
                    c+=1
                t+=1
            p+=1
            
        # Populate zero matrix
        zero_matrix = copy.deepcopy(new_X)
        p = 0
        for ptcp in stats_all_ptcps[stat_to_use].keys(): # These do not need to be arranged by time/space, it's just participant by participant
            t = 0
            for time in sorted(stats_all_ptcps[stat_to_use][ptcp]): # These are temporally arranged
                c = 0
                for chan in stats_all_ptcps[stat_to_use][ptcp][time].keys(): # These are not spatially arranged
                    zero_matrix[p][t][c] = 0.0
                    c+=1
                t+=1
            p+=1


     
        #### MAIN CHOICE AS TO WHAT TO PLOT
        #headshape = copy.deepcopy(self.group_A_avg['B'])  
        self.get_child_adult_sample()

        headshape = sample_individual.epochs_ransac_autoreject[sample_individual.cond_B] # (self.group_A_avg['B'])
        # headshape = centre_sensor_locations(headshape)[0]
        info = headshape.info  # self.group_A_avg['B'].info # diff_evokeds_group.info   # all_evoked_condition_A[0].info
        pos = mne.find_layout(info).pos
        pos_l = list(pos)
        # Delete location
        chan_nums_to_delete = sorted([int(x.replace("MEG ","")) for x in self.missing_channels],reverse=True)
        print("Chans to delete: ", chan_nums_to_delete)        
        #pos_l[(chan_nums_to_delete[0]-1):(chan_nums_to_delete[0]+1)]
        for index in chan_nums_to_delete:
            del pos_l[index]
        pos_l = pos_l[0:NUM_CHANS_TO_KEEP]        
        # After
        pos = np.array(pos_l)
        print("Number of channels not bad :", len(pos))



        # configure variables for visualization
        colors = {"rel_target": "crimson", "unrel_target": 'steelblue'}
        # get sensor positions via layout
        # find offset centroid and remove
        x = [p[0] for p in pos[:,0:2]]
        y = [p[1] for p in pos[:,0:2]]
        centroid = (sum(x) / len(pos), sum(y) / len(pos))
        # Scale by a factor of 5 (this is weird)
        pos_scaled_centred = np.array([[(points[0]-centroid[0])/5,(points[1]-centroid[1])/5,points[2],points[3]] for points in pos])
        
 
        # Adjacency matrix
        ch_adjacency, ch_names = mne.channels.find_ch_adjacency(info, 'mag');
        
        t_clust, clusters, p_values, H0 =  mne.stats.spatio_temporal_cluster_test(
            np.array([new_X,zero_matrix]),
            n_jobs=1,
            #threshold=threshold,
            adjacency=ch_adjacency,
            n_permutations=N_PERMUTATIONS
            ) 


        print("P values...", p_values)
        good_cluster_inds = np.where(p_values < CLUSTER_CUTOFF)[0]
        print("Good cluster indices: ", good_cluster_inds)
        
        filename_all_clusters = 'SpatioTemporalCluster '+stat_to_use + " ERFtoPredictor, " + waveform + " #All Details"
        filename_all_clusters = filename_all_clusters+ " "+timestring+file_name_addition+'.txt'
        with open(filename_all_clusters, 'w') as g:
            g.write("Good cluster indices\n")
            g.write(str(good_cluster_inds))
            g.write("\n")
            g.write("All cluster p values\n")
            g.write(str(p_values))
            g.write("\n")
            g.write("\n")            
            g.write("t-values of clusters\n")
            g.write(str(t_clust))   

        times = statistics_head.to_data_frame()['time'].values
        for i_clu, clu_idx in enumerate(good_cluster_inds):
            # unpack cluster information, get unique indices
            time_inds, space_inds = np.squeeze(clusters[clu_idx])

            ch_inds                 = np.unique(space_inds)
            time_inds               = np.unique(time_inds) 
            significant_channels    = ["MEG "+str(ch_ind).zfill(3) for ch_ind in ch_inds]
            significant_times       = [times[r] for r in time_inds]
            print("Significant channels ", str(significant_channels))
            print("Significant times", str(significant_times) ) #[sample_time+5*time_ind for time_ind in time_inds])

            # get topography for F stat
            t_map = t_clust[time_inds, ...].mean(axis=0)

            # get topography of difference
            time_shift = self.participants[0].epochs_ransac_autoreject.time_as_index(self.participants[0].epochs_ransac_autoreject.times[0])      # fix windowing shift
            diff_topo = np.mean(statistics_head.data[:,time_inds+time_shift],axis=1)   # statistics_head.data # 

            # find the time points of significance
            sig_times = self.participants[0].epochs_ransac_autoreject.times[time_inds]

            # create spatial mask
            mask = np.zeros((t_map.shape[0], 1), dtype=bool)
            mask[ch_inds, :] = True

            # initialize figure
            fig, ax_topo = plt.subplots(nrows=1, ncols=1, figsize=(20, 6)) # ncols=2,
            fig,(ax1,ax2) = plt.subplots(ncols=2)

            # plot average difference and mark significant sensors
            image, _ = plot_topomap(diff_topo,
                                    pos_scaled_centred,
                                    mask=mask, axes=ax_topo, cmap='RdBu_r', 
                                    # title=None,
                                    # sphere=(shift_x, shift_y, shift_z, CONSTANT_RADIUS_MULT*find_radius(statistics_head)),  # ADDED
                                    res=RESOLUTION,
                                    vmin=-1*10**(-14), vmax=1*10**(-14),
                                    #vmin=np.min, vmax=np.max
                                    show=False)

            # create additional axes (for ERF and colorbar)
            divider = make_axes_locatable(ax_topo)

            # add axes for colorbar
            ax_colorbar = divider.append_axes('right', size='10%', pad=0.05)
            plt.colorbar(image, cax=ax_colorbar)
            xlabel = "Relationship between ERF and predictor at time  "

            ax_topo.set_xlabel(xlabel+'({:0.3f} - {:0.3f} s)'.format(*sig_times[[0, -1]]), fontsize=FONTSIZE_LABELS)
            if stat_to_use == 'corr':
                ax_topo.set_ylabel("Correlation", fontsize=FONTSIZE_LABELS)
            else:
                ax_topo.set_ylabel(stat_to_use, fontsize=FONTSIZE_LABELS)
            
            # add new axis for time courses and plot time courses
            ax_signals = divider.append_axes('right', size='300%', pad=1.2)
            #ax_signals.set_xlim([times[0]*1/SEC_TO_MS, times[-1]*1/SEC_TO_MS+EPSILON_TIME])
            #ax_signals.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS, xmax=EPOCH_END_ANALYSIS*SEC_TO_MS)
            ax_signals.set_xlim(xmin=PLOT_START*SEC_TO_MS, xmax=PLOT_END*SEC_TO_MS)


            #print("Ax signals ", type(ax_signals), ax_signals)
            title = 'Cluster #{0}, {1} sensor'.format(i_clu + 1, len(ch_inds))
            if len(ch_inds) > 1:
                title+="s"
            title += " had a " + stat_to_use + " significantly different to zero"


            figure = plot_compare_evokeds([statistics_head],
                                          title=title,
                                          picks=ch_inds,
                                          axes=None, # ax_signals,
                                          #colors=colors,
                                          show= True, #False,
                                          ci=True,
                                          #cmap = None,
                                          show_sensors=True,
                                          legend=True,
                                          combine = 'mean',
                                          #split_legend=True,
                                          #truncate_yaxis='max_ticks'
                                         )

            #figure(figsize=(20, 6), dpi=240)
            # plot temporal cluster extent
            ymin, ymax = ax_signals.get_ylim()
            ax_signals.fill_betweenx((ymin, ymax), sig_times[0], sig_times[-1],
                                    color='orange', alpha=0.3)

            # clean up viz
            mne.viz.tight_layout(fig=fig)
            fig.subplots_adjust(bottom=.05)
            plt.legend('',frameon=False)      

            filename = 'SpatioTemporalCluster '+stat_to_use + " ERFtoPredictor, " + waveform + ' #'+str(clu_idx+1)
            figure[0].savefig(filename+ " "+timestring+file_name_addition+'.jpg',
                        format='jpeg',
                        dpi=DPI,
                        bbox_inches='tight')    
            with open(filename+" "+timestring+file_name_addition+' Details.txt', 'w') as f:
                f.write("Significant channels\n")
                f.write(str(significant_channels))
                f.write("\n")                     
                f.write("Significant times\n")
                f.write(str(significant_times))
                f.write("\n")                
                # f.write("t-map\n")
                # f.write(str(t_map))
                # f.write("\n")                                            
                # f.write("Mask\n")
                # f.write(str(mask))


        return return_val, n





    def adjust_adult_head_info(self):

        ## Get child head sample object
        self.get_child_adult_sample()
        self.sample_child = centre_sensor_locations(self.sample_child)[0]

        condition_A = self.sample_child.cond_A
        condition_B = self.sample_child.cond_B
        condition_names = [condition_A,condition_B]    

        
        # Replace all adult co-ordinates with those of children at mapped sensors
        for ptcp in self.participants:
            if ptcp.is_adult_system: 

                mne.channels.find_layout(ptcp.epochs_ransac_autoreject.info, ch_type='mag', exclude='bads')
                # mne.viz.plot_sensors(ptcp.epochs_ransac_autoreject.info,
                #                      kind='topomap',
                #                     #show_names = True
                #                     );

                # mne.viz.plot_evoked_topomap(ptcp.epochs_ransac_autoreject.average(),
                #                         #sphere=(shift_x, shift_y, shift_z, find_radius(self.grand_average_cond_A)),#,radius=radius_youngest)),
                #                         res=RESOLUTION, size=SIZE, cbar_fmt='%5.3f',
                #                         title = 'ERFs Before Adjustment',
                #                         #vmax=YLIM
                #                         outlines='head' # 'skirt'
                #                        );
  

                # Plot after centring, before mapping
                # mne.viz.plot_sensors(ptcp.epochs_ransac_autoreject.info,
                #                      kind='topomap',
                #                     #show_names = True
                #                     );        
                # avg = ptcp.epochs_ransac_autoreject.average()
                # shifts = centre_sensor_locations(avg)
                # shift_x,shift_y,shift_z = shifts[1], shifts[2], 0
                # mne.viz.plot_evoked_topomap(avg,
                #                     sphere=(shift_x, shift_y, shift_z, find_radius(avg)),#,radius=radius_youngest)),
                #                     res=RESOLUTION, size=SIZE, cbar_fmt='%5.3f',
                #                     title = 'Before mapping, after centring',
                #                     #vmax=YLIM
                #                     outlines='head'
                #                    );



                # Mapping
                #for adult_channel in ptcp.sensor_mapping_dct_simple.keys():
                ptcp.epochs_ransac_autoreject.info = self.sample_child.epochs_ransac_autoreject.info
                # m = 0
                # for adult_equivalent_channel in ptcp.epochs_ransac_autoreject.info['ch_names']:
                #     print("Adult channel ", adult_equivalent_channel)
                #     for ch in ptcp.sensor_mapping_dct_simple.keys():
                #         if ptcp.sensor_mapping_dct_simple[ch] == adult_equivalent_channel:
                #             child_equivalent_channel = ch
                #     print("Child equivalent channel ", child_equivalent_channel)
                #     p = 0
                #     for ch in self.sample_child.epochs_ransac_autoreject.info['chs']:
                #         if ch['ch_name'] == child_equivalent_channel: # self.sample_child.epochs_ransac_autoreject.info['chs'][p]
                #             #equiv_loc = ch['ch_name']
                #             #print("Equiv loc for ", adult_channel, " is ch ", equiv_loc )
                #             new_x, new_y = self.sample_child.epochs_ransac_autoreject.info['chs'][p]['loc'][0], self.sample_child.epochs_ransac_autoreject.info['chs'][p]['loc'][1]
                #             print(adult_equivalent_channel, child_equivalent_channel, m, p, new_x, new_y)

                #             ptcp.epochs_ransac_autoreject.info['chs'][m]['loc'][0] = new_x
                #             ptcp.epochs_ransac_autoreject.info['chs'][m]['loc'][1] = new_y
                #         p+=1
                #     m+=1




                # Plot before recentre 
                # mne.viz.plot_sensors(ptcp.epochs_ransac_autoreject.info,
                #                      kind='topomap',
                #                     show_names = True
                #                     );     
                # avg = ptcp.epochs_ransac_autoreject.average()
                # shifts = centre_sensor_locations(avg)
                # shift_x,shift_y,shift_z = shifts[1], shifts[2], 0    
                # mne.viz.plot_evoked_topomap(avg,
                #                     res=RESOLUTION, size=SIZE, cbar_fmt='%5.3f',
                #                     sphere=(shift_x, shift_y, shift_z, find_radius(avg)),#,radius=radius_youngest)),
                #                     title = 'After mapping, before re-centre',
                #                     #vmax=YLIM
                #                     outlines='head'
                #                    );




                # Re-centring
                ptcp.epochs_ransac_autoreject = centre_sensor_locations(ptcp.epochs_ransac_autoreject)[0]
                # Plot after recentre
                # avg = ptcp.epochs_ransac_autoreject.average()
                # shifts = centre_sensor_locations(avg)
                # shift_x,shift_y,shift_z = shifts[1], shifts[2], 0
                # mne.viz.plot_evoked_topomap(avg,
                #                     sphere=(shift_x, shift_y, shift_z, find_radius(avg)),#,radius=radius_youngest)),
                #                     res=RESOLUTION, size=SIZE, cbar_fmt='%5.3f',
                #                     title = 'After mapping and re-centre',
                #                     #vmax=YLIM
                #                     outlines='head'
                #                    );
                # mne.viz.plot_sensors(ptcp.epochs_ransac_autoreject.info,
                #                      kind='topomap',
                #                     show_names = True
                #                     );






    def all_stats_loop(self,tgt_dir,MMR_Experiment=None,age_group="Both",stat_to_use='corr',percentile_cutoff=20):
        
        '''
        For each experiment
            For unstandardised and standardised separately
                Generates GFP and ERF graphs and times of significant differences
                Generates spatiotemporal ERF difference graphs showing times of significant differences
                Generates statistical heads (e.g. for correlation)
                Runs regressions 

        stat_to_use choose from ['corr', 'slope']
        age_group choose from ['Both' (default), 'Young', 'Old']
        tgt_dir is the directory to source .pickles from (and save stats results to)

        percentile is for comparing low vs high surprise (and young vs old)
        '''
        


        problem_folders = []
        error_messages  = []
        folder_num      = 1
        if not os.path.isfile(tgt_dir):
            # Set up experiment
            # Add participants
            print("Target dir, ", tgt_dir, " folder #: ", folder_num)
            os.chdir(tgt_dir)
            file_names = glob.glob(f"{tgt_dir}*.pickle")
            if len(file_names) == 0:
                file_names = glob.glob(f"{tgt_dir}*Bins.pickle")

            # Figure out which experiment it is for
            with open(file_names[0], "rb") as f:
                ptcp = pickle.load(f)
            condition_to_compare = ptcp.cond_B
            events_to_tag = ptcp.events_to_tag  
            print("Events to tag ", events_to_tag)
            print("Condition B ", conds_to_compare[events_to_tag[1]][1])
            print("Condition to compare ", condition_to_compare)

            gc.collect()        
            for file in file_names:
                include = False
                for child_string in self.child_participant_strings:
                    if child_string in file:
                        include = True
                for adult_string in self.adult_participant_strings:
                    if adult_string in file:
                        include = True            
                if include:
                    wait_until_memory_free(required_memory = 3, max_wait_time_mins = 1.5) # This requires some memory (set conservatively here as I may run many threads, should only be 500MB)        
                    self.add_participants_from_disk([file])
                    ptcp = self.participants[-1]
                    if ptcp.cond_A in ptcp.evoked_generic.keys() and ptcp.cond_B in ptcp.evoked_generic.keys():
                        self.condition_to_compare = ptcp.cond_B
                        ptcp.epochs_ransac_autoreject.crop(tmin=EPOCH_START_ANALYSIS, tmax=EPOCH_END_ANALYSIS, include_tmax=True)
                        ptcp.evoked_generic[ptcp.cond_A].crop(tmin=EPOCH_START_ANALYSIS, tmax=EPOCH_END_ANALYSIS, include_tmax=True)
                        ptcp.evoked_generic[ptcp.cond_B].crop(tmin=EPOCH_START_ANALYSIS, tmax=EPOCH_END_ANALYSIS, include_tmax=True)
                        ptcp.evoked_all.crop(tmin=EPOCH_START_ANALYSIS, tmax=EPOCH_END_ANALYSIS, include_tmax=True)
                        del_attributes(ptcp, attrs_to_delete)
                        self.participants[-1] = copy.deepcopy(ptcp)                          
                    else:
                        self.warn("CONDITIONS ARE WRONG" +str(folder) +str(ptcp.p_id))
                        del self.participants[-1]
            
            self.removed = {}
            self.remove_unwanted_data() 
            #self.remove_bad_sound_delay()

            #self.align_channels()

            #self.adjust_adult_head_info()         
            
            age_num_groupings = int(100/20)
            self.find_age_groupings(num_age_groupings=age_num_groupings)
            self.group_participants_on_age(num_per_group=None, absolute=False, system=None, age_cutoff=20)
            self.compare_group(condition='B') 
            
    
            timestring = str(int(self.participants[0].epochs_ransac_autoreject.times[0]*SEC_TO_MS+EPSILON_TIME))+ "-"+str(int(self.participants[0].epochs_ransac_autoreject.times[-1]*SEC_TO_MS+EPSILON_TIME))+"ms"
            
            #self.compare_ERFs_conditions()
            self.cond_A_evoked = []
            self.cond_B_evoked = []
            r = 0
            for ptcp in self.participants:
                if r % 10 == 0:
                    print("Compare ERFs loop: participant #: ", r)
                c_a = ptcp.evoked_generic[ptcp.cond_A]
                c_b = ptcp.evoked_generic[ptcp.cond_B]
                # Somewhat meaningless 'sum' column
                c_a_df = sum_df(c_a.to_data_frame())
                c_b_df = sum_df(c_b.to_data_frame())
                sum_a = c_a_df['sum'].values
                sum_b = c_b_df['sum'].values
                if sum([1 for x in sum_a if math.isnan(x)]) >0 or sum([1 for x in sum_b if math.isnan(x)]) >0:
                    self.warn(" nans present for " + str(self.p_id) + ", not adding @@@@@@@@@@@@@@ ")
                else:
                    self.cond_A_evoked.append(c_a) 
                    self.cond_B_evoked.append(c_b)
                for var_to_delete in ['c_a','c_b']: #",'x','y']:
                    globals().pop(var_to_delete, None)
                gc.collect()
                r+=1
            self.grand_average_cond_A = mne.grand_average(self.cond_A_evoked)
            self.grand_average_cond_B = mne.grand_average(self.cond_B_evoked)
            self.grand_average_cond_A_df = self.grand_average_cond_A.to_data_frame()
            self.grand_average_cond_B_df = self.grand_average_cond_B.to_data_frame()
            #self.grand_average_cond_AB_diff = mne.combine_evoked([self.grand_average_cond_B,self.grand_average_cond_A],weights=[1,-1])
            self.grand_average_cond_AB_diff = self.grand_average_cond_B_df - self.grand_average_cond_A_df
            self.grand_average_cond_AB_diff['time'] = self.grand_average_cond_A_df['time'].values
            times = self.grand_average_cond_A_df['time'].values


        
            #Plot after recentre
            #shifts = centre_sensor_locations(self.grand_average_cond_A)
            #shift_x,shift_y,shift_z = shifts[1], shifts[2], 0
            mne.viz.plot_evoked_topomap(self.grand_average_cond_A,
                                #sphere=(shift_x, shift_y, shift_z, find_radius(self.grand_average_cond_A)),#,radius=radius_youngest)),
                                res=RESOLUTION, size=SIZE, cbar_fmt='%5.3f',
                                title = 'Condition A',
                                #vmax=YLIM
                                outlines='head'
                               );
            #shifts = centre_sensor_locations(self.grand_average_cond_B)
            #shift_x,shift_y,shift_z = shifts[1], shifts[2], 0
            mne.viz.plot_evoked_topomap(self.grand_average_cond_B,
                                #sphere=(shift_x, shift_y, shift_z, find_radius(self.grand_average_cond_B)),#,radius=radius_youngest)),
                                res=RESOLUTION, size=SIZE, cbar_fmt='%5.3f',
                                title = 'Condition B',
                                #vmax=YLIM
                                outlines='head'
                               );
        
        

            #try:        
            os.chdir(tgt_dir)


            ##### ERFs
            # Export df     time X channels
            #self.grand_average_cond_A_df = self.grand_average_cond_A_df.reindex(sorted(df.columns), axis=1)

            filename = "ERF Cond A "+str(percentile_cutoff)+"%ile "+timestring+".csv"
            self.grand_average_cond_A_df.to_csv(filename,index=False)     
            filename = "ERF Cond B "+str(percentile_cutoff)+"%ile "+timestring+".csv"
            self.grand_average_cond_B_df.to_csv(filename,index=False)          
            filename = "ERF Cond B minus Cond A "+str(percentile_cutoff)+"%ile "+timestring+".csv"
            self.grand_average_cond_AB_diff.to_csv(filename,index=False)      




            #### GFPs
            # Compute GFPs for each condition (using group evoked response)
            # And compute differences in GFP between each condition
            gfp_cond_A = RMS_DF(self.grand_average_cond_A_df)
            gfp_cond_A_df = pd.DataFrame({'time': self.grand_average_cond_A_df['time'].values, 'GFP_COND_A': gfp_cond_A})
            gfp_cond_B = RMS_DF(self.grand_average_cond_B_df)
            gfp_cond_B_df = pd.DataFrame({'time': self.grand_average_cond_B_df['time'].values, 'GFP_COND_B': gfp_cond_B})
            gfp_diff_conditions =  [gfp_cond_B[x] - gfp_cond_A[x] for x in range(0,len(gfp_cond_B))]
            times = np.arange(EPOCH_START_ANALYSIS*SEC_TO_MS,1+EPOCH_END_ANALYSIS*SEC_TO_MS,int((EPOCH_END_ANALYSIS*SEC_TO_MS-EPOCH_START_ANALYSIS*SEC_TO_MS)/(len(gfp_cond_A)-1)))
            #gfp_diff_conditions_df = pd.DataFrame([self.grand_average_cond_B_df['time'].values,gfp_diff_conditions],columns=['time','GFP_DIFF']  
            gfp_diff_conditions_df = pd.DataFrame({'time': self.grand_average_cond_A_df['time'].values, 'GFP_COND_(B-A)': gfp_diff_conditions})
            



            plt = plot_module(  xs = [times],
                                ys = [gfp_diff_conditions],
                                xlabel = "Time (ms)", ylabel="GFP Cond B - Cond A",
                                legends=["GFP Cond B - Cond A"],
                                title="High surprise RMS\n"+"Minus low surprise RMS",
                                font_object = self.font,
                                x_min_to_show = PLOT_START*SEC_TO_MS,
                                x_max_to_show = PLOT_END*SEC_TO_MS,                                  
                                #x_min_to_show = EPOCH_START_ANALYSIS*SEC_TO_MS,
                                #x_max_to_show = EPOCH_END_ANALYSIS*SEC_TO_MS                            
                                )
            # plt.plot(times,gfp_diff_conditions)
            # #plt.xlim([min(times),max(times)])
            # plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed 
            # plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
            # plt.title("High surprise RMS minus low surprise RMS", fontsize=FONTSIZE_TITLE)
            plt.show()



            filename = "GFP Cond A "+str(percentile_cutoff)+"%ile "+timestring+".csv"
            gfp_cond_A_df.to_csv(filename,index=False)      
            filename = "GFP Cond B "+str(percentile_cutoff)+"%ile "+timestring+".csv"
            gfp_cond_B_df.to_csv(filename,index=False)      
            filename = "GFP Cond B minus Cond A "+str(percentile_cutoff)+"%ile "+timestring+".csv"
            gfp_diff_conditions_df.to_csv(filename,index=False)  
  


            # Construct DF with all the GFP comparisons
            df_gfps_conditions = pd.DataFrame({'time': self.grand_average_cond_A_df['time'].values,
                                    'GFP CondA': gfp_cond_A,
                                    'GFP CondB': gfp_cond_B,
                                    'GFP CondB-A': gfp_diff_conditions,
                                               }
                                              )
            df_gfps_conditions.to_csv("GFP Conditions "+str(percentile_cutoff)+"%ile "+timestring+".csv",index=False)      
            # Statistically test the difference between conditions using each individual's GFP
            self.save_gfp_diff_cond(filter_age_max=6)




            # Plot differences together with raw values
            fig,ax = plt.subplots()
            # make a plot
            NUM_TICKS = 5
            y_min = min(min(df_gfps_conditions['GFP CondA'].values), min(df_gfps_conditions['GFP CondB'].values))
            y_max = max(max(df_gfps_conditions['GFP CondA'].values), max(df_gfps_conditions['GFP CondB'].values))
            ax.set_ylim(ymin=y_min,ymax=y_max) 
            y_ticks = np.array([y_min+(y_max-y_min)*(r)/(NUM_TICKS-1) for r in range(0,NUM_TICKS)])   
            ax.set_yticks(y_ticks)
            #ax.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS, xmax=EPOCH_END_ANALYSIS*SEC_TO_MS)   
            ax.set_xlim(xmin=PLOT_START*SEC_TO_MS,xmax=PLOT_END*SEC_TO_MS)


            ax.plot(df_gfps_conditions['time'].values, df_gfps_conditions['GFP CondA'].values, label='Cond A', color="red") # ,marker="o"
            ax.plot(df_gfps_conditions['time'].values, df_gfps_conditions['GFP CondB'].values, label='Cond B', color="blue" )
            # Hide the right and top spines
            #ax.spines['right'].set_visible(False)
            #ax.spines['top'].set_visible(False)
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')
            ax.legend(bbox_to_anchor=(1, 0.5), fontsize=FONTSIZE_LEGEND)
            # set x-axis label
            ax.set_xlabel("Time (ms)",   fontsize=FONTSIZE_LABELS)
            # set y-axis label
            ax.set_ylabel("GFP (T)",     fontsize=FONTSIZE_LABELS) # color="red",
            fig.savefig('GFP vs condition One Only '+timestring+'.jpg',
                        format='jpeg',
                        dpi=DPI,
                        bbox_inches='tight')

            # twin object for two different y-axis on the sample plot
            ax2=ax.twinx()
            # make a plot with different y-axis using second axis object
            #ax2.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS, xmax=EPOCH_END_ANALYSIS*SEC_TO_MS)   
            ax2.set_xlim(xmin=PLOT_START*SEC_TO_MS,xmax=PLOT_END*SEC_TO_MS)
            ax2.plot(df_gfps_conditions['time'].values, df_gfps_conditions['GFP CondB-A'].values, label='B-A', color="green")
            # Hide the right and top spines
            #ax2.spines.right.set_visible(False)
            #ax2.spines.top.set_visible(False)
            # Only show ticks on the left and bottom spines
            #ax2.yaxis.set_ticks_position('left')
            #ax2.xaxis.set_ticks_position('bottom')
            ax2.legend(fontsize=FONTSIZE_LEGEND, bbox_to_anchor=(1, 0.5) )
            ax2.set_ylabel("GFP of [Grand avg Condition B minus A]",color="green",fontsize=FONTSIZE_LABELS)

            #plt.legend([("GFP Condition A "), ("GFP Condition B "), ("GFP Condition B minus Condition A")], fontsize=FONTSIZE_LEGEND)
            plt.title("GFP for different conditions", fontsize=FONTSIZE_TITLE)
            # save the plot as a file
            fig.savefig('GFP vs condition '+timestring+'.jpg',
                        format='jpeg',
                        dpi=DPI,
                        bbox_inches='tight')
            plt.clf()
            # Plot condition differences only



            plt = plot_module(  xs = [df_gfps_conditions['time'].values],
                                ys = [[df_gfps_conditions['GFP CondB'].values[x] - df_gfps_conditions['GFP CondA'].values[x] for x in range(0,len(df_gfps_conditions['GFP CondA'].values))]],
                                xlabel = "Time (ms)", ylabel="GFP difference (T)",
                                legends=["GFP Cond B\n"+"minus GFP Cond A"],
                                title="GFP of (Grand avg) Condition B minus\n"+"GFP of (Grand avg) Condition A",
                                font_object = self.font,
                                colours = [], # 
                                x_min_to_show = PLOT_START*SEC_TO_MS,
                                x_max_to_show = PLOT_END*SEC_TO_MS,                                  
                                #x_min_to_show = EPOCH_START_ANALYSIS*SEC_TO_MS,
                                #x_max_to_show = EPOCH_END_ANALYSIS*SEC_TO_MS                            
                                )

            # fig=plt.figure()
            # ax=fig.add_subplot(111)
            # ax.plot(df_gfps_conditions['time'].values, [df_gfps_conditions['GFP CondB'].values[x] - df_gfps_conditions['GFP CondA'].values[x] for x in range(0,len(df_gfps_conditions['GFP CondA'].values))], label = "GFP Cond B\n"+str("minus GFP Cond A"))  # plt.
            # ax.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS,xmax=EPOCH_END_ANALYSIS*SEC_TO_MS)      # Makes first x tick start at origin 
            # ax.set_xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)  
            # ax.set_ylabel("GFP difference (T)", fontsize=FONTSIZE_LABELS, labelpad=PADDING, fontname="Times New Roman")    
            # ax.legend(loc='best', fontsize=FONTSIZE_LEGEND, bbox_to_anchor=(1, 0.5), prop=self.font )

            # #plt.plot(df_gfps_conditions['time'].values, [df_gfps_conditions['GFP CondB'].values[x] - df_gfps_conditions['GFP CondA'].values[x] for x in range(0,len(df_gfps_conditions['GFP CondA'].values))])
            # #plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed 
            # plt.title("GFP of (Grand avg of) Condition B minus A", fontsize=FONTSIZE_TITLE)
            # #plt.ylabel("GFP difference (T)", fontsize=FONTSIZE_LABELS)
            # #plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
            # #plt.legend(["GFP Cond B\n"+str("minus GFP Cond A")], fontsize=FONTSIZE_LEGEND, loc='best', bbox_to_anchor=(1, 0.5))

            plt.savefig('GFP vs condition, differences '+timestring+'.jpg', dpi=DPI)
            plt.clf()






            timestring = str(int(self.participants[0].epochs_ransac_autoreject.times[0]*SEC_TO_MS+EPSILON_TIME))+ "-"+str(int(self.participants[0].epochs_ransac_autoreject.times[-1]*SEC_TO_MS+EPSILON_TIME))+"ms"

            #### @@@@@ BEGIN AGE ANALYSIS
            # Youngest vs oldest 10%
            self.group_analysis(size_age_bucket_percent=10);
            self.save_gfp_diff_ages(size_age_bucket_percent=10)
            df_gfps_group_10 = self.export_group(size_age_bucket_percent=10)
            df_gfps_group_10.to_csv("GFP Group 10"+" %ile "+timestring+".csv",index=False)    

            # Youngest vs oldest 20%
            self.group_analysis(size_age_bucket_percent=20);   
            self.save_gfp_diff_ages(size_age_bucket_percent=20)        
            df_gfps_group_20 = self.export_group(size_age_bucket_percent=20)
            df_gfps_group_20.to_csv("GFP Group 20"+" %ile "+timestring+".csv",index=False)    
            #shifts = centre_sensor_locations(self.group_A_avg['B'])
            #shift_x,shift_y,shift_z = shifts[1], shifts[2], 0
            mne.viz.plot_evoked_topomap(self.group_A_avg['B'],
                                #sphere=(shift_x, shift_y, shift_z, find_radius(self.group_A_avg['B'])),#,radius=radius_youngest)),
                                res=RESOLUTION, size=SIZE, cbar_fmt='%5.3f',
                                title = 'Group A Condition B',
                                #vmax=YLIM
                                outlines='head'
                               );
            #shifts = centre_sensor_locations(self.group_B_avg['B'])
            #shift_x,shift_y,shift_z = shifts[1], shifts[2], 0
            mne.viz.plot_evoked_topomap(self.group_B_avg['B'],
                                #sphere=(shift_x, shift_y, shift_z, find_radius(self.group_B_avg['B'])),#,radius=radius_youngest)),
                                res=RESOLUTION, size=SIZE, cbar_fmt='%5.3f',
                                title = 'Group B Condition B',
                                #vmax=YLIM
                                outlines='head'
                               );                                  
                                                                        
                                      
            # Use specific age groupings
            self.group_analysis(age_bounds_low=age_bounds_low,age_bounds_high=age_bounds_high) # size_age_bucket_percent=20,
            self.save_gfp_diff_ages(size_age_bucket_percent=20,override_filename=str(age_bounds_low)+str(age_bounds_high))
            df_gfps_group_custom = self.export_group(size_age_bucket_percent=20)
            str_ = str(age_bounds_low)+" and "+str(age_bounds_high)
            df_gfps_group_custom.to_csv("GFP Group ages "+str_+' '+timestring+'.csv',index=False)      

            ## Plot oldest vs youngest together
            # Method A

            NUM_TICKS = 3
            x_ticks = df_gfps_group_10['time'].values
            y_min = min(min(df_gfps_group_10['GFP Group A EveryTrial'].values),
                        min(df_gfps_group_20['GFP Group A EveryTrial'].values),
                        min(df_gfps_group_10['GFP Group B EveryTrial'].values),
                        min(df_gfps_group_20['GFP Group B EveryTrial'].values)
                        )
            y_max = max(max(df_gfps_group_10['GFP Group A EveryTrial'].values),
                        max(df_gfps_group_20['GFP Group B EveryTrial'].values),
                        max(df_gfps_group_10['GFP Group B EveryTrial'].values),
                        max(df_gfps_group_20['GFP Group B EveryTrial'].values)
                        )   

            y_ticks = np.array([y_min+(y_max-y_min)*(r)/(NUM_TICKS-1) for r in range(0,NUM_TICKS)])
            # plt = plot_module(  xs = [x_ticks,x_ticks,x_ticks,x_ticks],
            #                     ys = [df_gfps_group_10['GFP Group A EveryTrial'].values,
            #                             df_gfps_group_10['GFP Group B EveryTrial'].values,
            #                             df_gfps_group_20['GFP Group A EveryTrial'].values,
            #                             df_gfps_group_20['GFP Group B EveryTrial'].values],
            #                     xlabel = "Time (ms)", ylabel="GFP (T)",
            #                     legends=["Youngest "+str(PERCENTILES_AGE[0])+ "%",
            #                             "Oldest "+str(PERCENTILES_AGE[0])+ " %",
            #                             "Youngest "+str(PERCENTILES_AGE[1])+ " %",
            #                             "Oldest "+str(PERCENTILES_AGE[1])+ "%"],
            #                     title="GFP for different age groups",
            #                     font_object = self.font,
            #                     xticklabels = x_ticks,
            #                     yticklabels = y_ticks,
            #                     colours = [], # 
            #                     x_min_to_show = PLOT_START*SEC_TO_MS,
            #                     x_max_to_show = PLOT_END*SEC_TO_MS,                                  
            #                     #x_min_to_show = EPOCH_START_ANALYSIS*SEC_TO_MS,
            #                     #x_max_to_show = EPOCH_END_ANALYSIS*SEC_TO_MS,
            #                     y_min_to_show = y_min,
            #                     y_max_to_show = y_max                                                             
            #                     )

    

            fig=plt.figure()
            ax=fig.add_subplot(111)
            ax.plot(x_ticks, df_gfps_group_10['GFP Group A EveryTrial'].values, 'b', linewidth=1, label="Youngest "+str(PERCENTILES_AGE[0])+ "%")
            ax.plot(x_ticks, df_gfps_group_10['GFP Group B EveryTrial'].values, 'b', linewidth=1, label="Oldest "+str(PERCENTILES_AGE[0])+ " %")
            ax.plot(x_ticks, df_gfps_group_20['GFP Group A EveryTrial'].values, 'b', linewidth=1, label="Youngest "+str(PERCENTILES_AGE[1])+ " %")
            ax.plot(x_ticks, df_gfps_group_20['GFP Group B EveryTrial'].values, 'b', linewidth=1, label="Oldest "+str(PERCENTILES_AGE[1])+ "%")        
            # Hide the right and top spines
            #ax.spines['top'].set_visible(False)
            #ax.spines['right'].set_visible(False)
            # Only show ticks on the left and bottom spines
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')
            ax.set_yticklabels(y_ticks,fontsize=FONTSIZE_AXES)
            plt.ylim([y_min,y_max])
            ax.set_yticks(y_ticks)        
            ax.set_ylabel("GFP (T)", fontsize=FONTSIZE_LABELS, labelpad=PADDING)
            ax.set_xticklabels(x_ticks, fontsize=FONTSIZE_AXES)
            #plt.xlim([min(df_gfps_group_10['time'].values),max(df_gfps_group_10['time'].values)])
            plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed 
            ax.set_xticks(x_ticks)            
            ax.set_xlabel('Time (ms)', fontsize=FONTSIZE_LABELS, labelpad=PADDING)
            ax.legend(loc='best', frameon=False, fontsize=FONTSIZE_LEGEND, bbox_to_anchor=(1, 0.5))
            # Method B
            # plt.plot(df_gfps_group_10['time'].values, df_gfps_group_10['GFP Group A EveryTrial'].values)
            # plt.plot(df_gfps_group_10['time'].values, df_gfps_group_10['GFP Group B EveryTrial'].values)
            # plt.plot(df_gfps_group_20['time'].values, df_gfps_group_20['GFP Group A EveryTrial'].values)
            # plt.plot(df_gfps_group_20['time'].values, df_gfps_group_20['GFP Group B EveryTrial'].values)
            #plt.legend(["Youngest "+str(PERCENTILES_AGE[0])+ "%", "Oldest "+str(PERCENTILES_AGE[0])+ " %", "Youngest "+str(PERCENTILES_AGE[1])+ " %", "Oldest "+str(PERCENTILES_AGE[1])+ "%",], fontsize=FONTSIZE_LEGEND)
            plt.title("GFP for different age groups", fontsize=FONTSIZE_TITLE)
            # plt.ylabel("GFP (T)", fontsize=FONTSIZE_LABELS)
            # plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
            plt.savefig('GFP vs age '+str(PERCENTILES_AGE[0])+'%-'+str(PERCENTILES_AGE[1])+"%"+timestring+'.jpg', dpi=DPI)
            # Image.open('testplot.png').save('testplot.jpg','JPEG')  
            plt.clf()




            # Plot oldest vs youngest differences
            plt = plot_module(  xs = [df_gfps_group_10['time'].values,df_gfps_group_20['time'].values],
                                ys = [
                                        [df_gfps_group_10['GFP Group A EveryTrial'].values[x] - df_gfps_group_10['GFP Group B EveryTrial'].values[x] for x in range(0,len(df_gfps_group_10['GFP Group B EveryTrial'].values))],
                                        [df_gfps_group_20['GFP Group A EveryTrial'].values[x] - df_gfps_group_20['GFP Group B EveryTrial'].values[x] for x in range(0,len(df_gfps_group_20['GFP Group B EveryTrial'].values))],
                                     ],
                                xlabel = "Time (ms)", ylabel="GFP difference (T)",
                                legends=["Youngest "+str(PERCENTILES_AGE[0])+ " % minus oldest "+str(PERCENTILES_AGE[0])+ "%",
                                        "Youngest "+str(PERCENTILES_AGE[1])+ " % minus oldest "+str(PERCENTILES_AGE[1])+ "%"
                                        ],
                                title="GFP difference between age groups",
                                font_object = self.font,
                                # xticklabels = x_ticks,
                                # yticklabels = y_ticks,
                                colours = [], # 
                                x_min_to_show = PLOT_START*SEC_TO_MS,
                                x_max_to_show = PLOT_END*SEC_TO_MS,                                 
                                # x_min_to_show = EPOCH_START_ANALYSIS*SEC_TO_MS,
                                # x_max_to_show = EPOCH_END_ANALYSIS*SEC_TO_MS,                                                           
                                )

            # fig=plt.figure()
            # ax=fig.add_subplot(111)
            # ax.plot(df_gfps_group_10['time'].values, [df_gfps_group_10['GFP Group A EveryTrial'].values[x] - df_gfps_group_10['GFP Group B EveryTrial'].values[x] for x in range(0,len(df_gfps_group_10['GFP Group B EveryTrial'].values))], label = "Youngest "+str(PERCENTILES_AGE[0])+ " % minus oldest "+str(PERCENTILES_AGE[0])+ "%")  # plt.
            # ax.plot(df_gfps_group_20['time'].values, [df_gfps_group_20['GFP Group A EveryTrial'].values[x] - df_gfps_group_20['GFP Group B EveryTrial'].values[x] for x in range(0,len(df_gfps_group_20['GFP Group B EveryTrial'].values))], label = "Youngest "+str(PERCENTILES_AGE[1])+ " % minus oldest "+str(PERCENTILES_AGE[1])+ "%")
            # ax.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS,xmax=EPOCH_END_ANALYSIS*SEC_TO_MS)      # Makes first x tick start at origin 
            # ax.set_xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)  
            # ax.set_ylabel("GFP difference (T)", fontsize=FONTSIZE_LABELS, labelpad=PADDING, fontname="Times New Roman")    
            # ax.legend(loc='best', fontsize=FONTSIZE_LEGEND, bbox_to_anchor=(1, 0.5), prop=self.font )

            # #plt.plot(df_gfps_group_10['time'].values, [df_gfps_group_10['GFP Group A EveryTrial'].values[x] - df_gfps_group_10['GFP Group B EveryTrial'].values[x] for x in range(0,len(df_gfps_group_10['GFP Group B EveryTrial'].values))])
            # #plt.plot(df_gfps_group_20['time'].values, [df_gfps_group_20['GFP Group A EveryTrial'].values[x] - df_gfps_group_20['GFP Group B EveryTrial'].values[x] for x in range(0,len(df_gfps_group_20['GFP Group B EveryTrial'].values))])
            # #plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed 
            # #plt.legend(["Youngest "+str(PERCENTILES_AGE[0])+ " % minus oldest "+str(PERCENTILES_AGE[0])+ "%",  "Youngest "+str(PERCENTILES_AGE[1])+ " % minus oldest "+str(PERCENTILES_AGE[1])+ "%"], loc='best', fontsize=FONTSIZE_LEGEND)
            # #plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
            # #plt.ylabel("GFP difference (T)",fontsize=FONTSIZE_LABELS)
            # plt.title("GFP of (Grand avg of) youngest minus oldest", fontsize=FONTSIZE_TITLE)
            plt.savefig('GFP vs age, differences '+timestring+'.jpg', dpi=DPI)
            # Image.open('testplot.png').save('testplot.jpg','JPEG')  
            plt.clf()


            # Plot custom age groups together
            plt = plot_module(  xs = [df_gfps_group_custom['time'].values,df_gfps_group_custom['time'].values],
                                ys = [
                                     df_gfps_group_custom['GFP Group A EveryTrial'].values,
                                     df_gfps_group_custom['GFP Group B EveryTrial'].values,
                                     ],
                                xlabel = "Time (ms)", ylabel="GFP (T)",
                                legends=["Age group " +str(age_bounds_low),
                                        "Age group " +str(age_bounds_high)
                                        ],
                                title="GFP for (Grand avg of)\n"+"two different age groups",
                                font_object = self.font,
                                # xticklabels = x_ticks,
                                # yticklabels = y_ticks,
                                colours = [], # 
                                x_min_to_show = PLOT_START*SEC_TO_MS,
                                x_max_to_show = PLOT_END*SEC_TO_MS,                                  
                                # x_min_to_show = EPOCH_START_ANALYSIS*SEC_TO_MS,
                                # x_max_to_show = EPOCH_END_ANALYSIS*SEC_TO_MS,                                                            
                                )

            # fig=plt.figure()
            # ax=fig.add_subplot(111)
            # ax.plot(df_gfps_group_custom['time'].values, df_gfps_group_custom['GFP Group A EveryTrial'].values, label = "Age group " +str(age_bounds_low))  # plt.
            # ax.plot(df_gfps_group_custom['time'].values, df_gfps_group_custom['GFP Group B EveryTrial'].values, label = "Age group " +str(age_bounds_high))
            # ax.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS,xmax=EPOCH_END_ANALYSIS*SEC_TO_MS)      # Makes first x tick start at origin 
            # ax.set_xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)  
            # ax.set_ylabel("GFP (T)", fontsize=FONTSIZE_LABELS, labelpad=PADDING, fontname="Times New Roman")    
            # ax.legend(loc='best', fontsize=FONTSIZE_LEGEND, bbox_to_anchor=(1, 0.5), prop=self.font )

            # #plt.plot(df_gfps_group_custom['time'].values, df_gfps_group_custom['GFP Group A EveryTrial'].values)
            # #plt.plot(df_gfps_group_custom['time'].values, df_gfps_group_custom['GFP Group B EveryTrial'].values)
            # #plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed 
            # #plt.legend(["Age group " +str(age_bounds_low), "Age group "+str(age_bounds_high) ], fontsize=FONTSIZE_LEGEND)
            # #plt.ylabel("GFP (T)", fontsize=FONTSIZE_LABELS)
            # #plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
            # plt.title("GFP for (Grand avg of)\n"+" two different age groups", fontsize=FONTSIZE_TITLE)            
            plt.savefig('GFP vs age '+timestring+'.jpg', dpi=DPI)
            # Image.open('testplot.png').save('testplot.jpg','JPEG')  
            plt.clf()


            # Plot custom age groups differences only
            # Plot custom age groups together


            plt = plot_module(  xs = [df_gfps_group_custom['time'].values],
                                ys = [
                                        [df_gfps_group_custom['GFP Group A EveryTrial'].values[x] - df_gfps_group_custom['GFP Group B EveryTrial'].values[x] for x in range(0,len(df_gfps_group_custom['GFP Group B EveryTrial'].values))],
                                     ],
                                xlabel = "Time (ms)", ylabel="GFP difference (T)",
                                legends=["Age group " +str(age_bounds_low)+"\n"+" minus age group "+str(age_bounds_high)
                                        ],
                                title="GFP for (Grand avg of)\n"+"two different age groups",
                                font_object = self.font,
                                # xticklabels = x_ticks,
                                # yticklabels = y_ticks,
                                colours = [], # 
                                x_min_to_show = PLOT_START*SEC_TO_MS,
                                x_max_to_show = PLOT_END*SEC_TO_MS,                                  
                                # x_min_to_show = EPOCH_START_ANALYSIS*SEC_TO_MS,
                                # x_max_to_show = EPOCH_END_ANALYSIS*SEC_TO_MS,                                                           
                                )

            # fig=plt.figure()
            # ax=fig.add_subplot(111)
            # ax.plot(df_gfps_group_custom['time'].values, [df_gfps_group_custom['GFP Group A EveryTrial'].values[x] - df_gfps_group_custom['GFP Group B EveryTrial'].values[x] for x in range(0,len(df_gfps_group_custom['GFP Group B EveryTrial'].values))], label = "Age group " +str(age_bounds_low) + " minus age group "+str(age_bounds_high))  # plt.
            # ax.set_xlim(xmin=EPOCH_START_ANALYSIS*SEC_TO_MS,xmax=EPOCH_END_ANALYSIS*SEC_TO_MS)      # Makes first x tick start at origin 
            # ax.set_xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)  
            # ax.set_ylabel("GFP difference (T)", fontsize=FONTSIZE_LABELS, labelpad=PADDING, fontname="Times New Roman")    
            # ax.legend(loc='best', fontsize=FONTSIZE_LEGEND, bbox_to_anchor=(1, 0.5), prop=self.font )

            # #plt.plot(df_gfps_group_custom['time'].values, [df_gfps_group_custom['GFP Group A EveryTrial'].values[x] - df_gfps_group_custom['GFP Group B EveryTrial'].values[x] for x in range(0,len(df_gfps_group_custom['GFP Group B EveryTrial'].values))])
            # #plt.xlim([EPOCH_START_ANALYSIS*SEC_TO_MS,EPOCH_END_ANALYSIS*SEC_TO_MS]) # Changed 
            # #plt.legend(["Age group " +str(age_bounds_low) + " minus age group "+str(age_bounds_high)], loc='best', fontsize=FONTSIZE_LEGEND)
            # #plt.ylabel("GFP difference (T)",fontsize=FONTSIZE_LABELS)
            # #plt.xlabel("Time (ms)", fontsize=FONTSIZE_LABELS)
            # plt.title("GFP for (Grand avg of)'n"+"two different age groups", fontsize=FONTSIZE_TITLE)            
            plt.savefig('GFP vs age '+str_+', differences '+timestring+'.jpg', dpi=DPI)
            # Image.open('testplot.png').save('testplot.jpg','JPEG')          
            plt.clf()



            self.cluster_analysis_spatio_temporal_v2(mode='condition',condition_comparison='B')                                         # Compare condition A vs condition B
            self.cluster_analysis_spatio_temporal_v2(mode='mmr_by_age',size_age_bucket_percent=10)                                      # Compare (B-A) vs age        
            self.cluster_analysis_spatio_temporal_v2(mode='mmr_by_age',size_age_bucket_percent=20)                                      # Compare (B-A) vs age
            self.cluster_analysis_spatio_temporal_v2(mode='age',condition_comparison=None, size_age_bucket_percent=10)                  # Compare all conditions across ages        
            self.cluster_analysis_spatio_temporal_v2(mode='age',condition_comparison=None, size_age_bucket_percent=20)                  # Compare all conditions across ages
            self.cluster_analysis_spatio_temporal_v2(mode='condition_by_age',condition_comparison='A', size_age_bucket_percent=10)      # Compare condition A across all ages
            self.cluster_analysis_spatio_temporal_v2(mode='condition_by_age',condition_comparison='B', size_age_bucket_percent=10)      # Compare condition B across all ages        
            self.cluster_analysis_spatio_temporal_v2(mode='condition_by_age',condition_comparison='A', size_age_bucket_percent=20)      # Compare condition A across all ages
            self.cluster_analysis_spatio_temporal_v2(mode='condition_by_age',condition_comparison='B', size_age_bucket_percent=20)      # Compare condition B across all ages
           



            method_dct = self.return_method_dct('On high surprise trials')
            assert method_dct['On high surprise trials']+method_dct['On low surprise trials']+method_dct['On every trial'] == 1, "Can only choose one option"
            waveform = self.return_waveform_str(method_dct)
            print("Waveform to analyse :", waveform)           
            corrs_high_exp_surprise_trials, n = self.statistical_head_spatio_temporal(method_dct,waveform,age_group=age_group,stat_to_use='corr') 
            slopes_high_exp_surprise_trials, n = self.statistical_head_spatio_temporal(method_dct,waveform,age_group=age_group,stat_to_use='slope')
            print("Number of high expected surprise trials ", n)

            method_dct = self.return_method_dct('On low surprise trials')
            assert method_dct['On high surprise trials']+method_dct['On low surprise trials']+method_dct['On every trial'] == 1, "Can only choose one option"
            waveform = self.return_waveform_str(method_dct)
            print("Waveform to analyse :", waveform)        

            corrs_low_exp_surprise_trials, n = self.statistical_head_spatio_temporal(method_dct,waveform,age_group=age_group,stat_to_use='corr') 
            slopes_low_exp_surprise_trials, n = self.statistical_head_spatio_temporal(method_dct,waveform,age_group=age_group,stat_to_use='slope') 
            print("Number of low expected surprise trials ", n)



            # Graphing
            filename = "Correlation predictor to GFP "+str(age_bounds_low)+"-"+str(age_bounds_high)+" "+str(age_group)
            corrs1 = corrs_low_exp_surprise_trials[0:]
            corrs2 = corrs_high_exp_surprise_trials[0:]
            threshold = abs(r_confidence(p_cutoff=0.05, n=n)[0])
            threshold_line = [threshold for x in corrs1]

            #print("H1",threshold_line)
            if (min(corrs1) < 0 and min(corrs2) <=0)  or ( (min(corrs1) < 0 and min(corrs2) >=0) or (min(corrs1) > 0 and min(corrs1) <=0)):
                threshold_line = [-x for x in threshold_line]   
                #filename="Thres Flipped "+filename
            filename = filename+".jpg"   
            
            corrs1f = ['{0:.3f}'.format(x) for x in corrs1]
            corrs2f = ['{0:.3f}'.format(x) for x in corrs2]
            threshold_linef = ['{0:.3f}'.format(x) for x in threshold_line]          
            corrs1f = [float(x) for x in corrs1f]
            corrs2f = [float(x) for x in corrs2f]
            threshold_linef = [float(x) for x in threshold_linef]

            padding_graph = 0.0 
            # y_min = min(min(corrs1f),min(corrs2f))
            # y_max = max(max(corrs1f),max(corrs2f))                       
            y_min = 0 if min(min(corrs1f),min(corrs2f)) >= 0 else (1+padding_graph)*min(min(corrs1f),min(corrs2f))
            y_max = 0 if max(max(corrs1f),max(corrs2f)) <= 0 else (1+padding_graph)*max(max(corrs1f),max(corrs2f))
            # y_min = y_min + (y_min-y_max)*padding_graph
            # y_max = y_max + (y_max-y_min)*padding_graph

            #times = [int(SEC_TO_MS*EPOCH_START_ANALYSIS + EPOCH_SIZE_MS*x) for x in range(0,(EPOCH_START_ANALYSIS-EPOCH_END_ANALYSIS)/EPOCH_SIZE_MS + 1)]

            #plt.rcParams["figure.figsize"] = (8,5.5)
            if stat_to_use == 'corr':
                ylabel='Correlation'
            elif stat_to_use == 'slope':
                ylabel = 'Slope'
            min_time = min(df_gfps_group_custom['time'].values)
            max_time = max(df_gfps_group_custom['time'].values)
            self.corr_low = corrs1f
            self.corr_high = corrs2f
            self.corr_thres = threshold_linef
            filename_new = filename.replace(".jpg","")
            filename_new = filename_new + " "+str(min_time)+"-"+str(max_time)+" ms.jpg" 
            compare_corrs_to_threshold(times, ylabel= ylabel, line1=corrs1f, line2=corrs2f, threshold_line=threshold_linef, y_min=y_min, y_max=y_max, filename=filename_new)
                
                  
                        
            folder_num+=1
    #         except Exception as e:
    #             print("PROBLEM ", folder, e)
    #             problem_folders.append(folder)
    #             error_messages.append(e)



        print("Problems with folders ", problem_folders)        
        

    def warn(self,warning):
        '''Print as a warning'''
        print("%%%%%%%%%%%%%%%% WARNING: "+ str(warning) + " %%%%%%%%%%%%%%%")

    def check_for_trigger_leakage(self):
        '''Used to check whether there is trigger leakage'''
        for ptcp in self.participants:
    
            ptcp.get_raw()
            ptcp.raw_df = ptcp.raw.to_data_frame()
            ptcp.stims = ptcp.raw_df[ptcp.raw_df["STI 014"] != 0.0]
            ptcp.stims['trigger_shift'] = ptcp.stims['STI 014'].shift(1)
            ptcp.stims['trigger'] = (ptcp.stims['trigger_shift']!=ptcp.stims['STI 014']) + 0
            ptcp.stims = ptcp.stims[ptcp.stims['trigger'] == 1]
            ptcp.stims['trigger_time'] = ptcp.stims['trigger']*ptcp.stims['time']
            
            if len(ptcp.stims) > 0:
                trigger_times = sorted(list(set(ptcp.stims['trigger_time'].values)))
                del ptcp.stims['trigger_shift']
                ptcp.stims = ptcp.stims.reset_index(drop=True)


                ptcp.epochs_surrounding = {}
                time_to_plot = 20 #ms
                r = 0
                for trigger_time in trigger_times:
                    for chan in ptcp.raw_df.columns.values:
                        if "MEG" in chan:
                            if chan in ptcp.epochs_surrounding.keys():
                                if trigger_time >= 1 and trigger_time < len(ptcp.raw_df)-time_to_plot:
                                    ptcp.epochs_surrounding[chan].append([ptcp.raw_df[chan][trigger_time-1]])    
                                    for time_delay in range(0,time_to_plot+1):
                                        ptcp.epochs_surrounding[chan][-1].append(ptcp.raw_df[chan][trigger_time+time_delay])
                            else:
                                if trigger_time >= 1 and trigger_time < len(ptcp.raw_df)-time_to_plot:
                                    ptcp.epochs_surrounding[chan] = []
                                    ptcp.epochs_surrounding[chan].append([ptcp.raw_df[chan][trigger_time-1]])  
                                    for time_delay in range(0,time_to_plot+1):
                                        ptcp.epochs_surrounding[chan][-1].append(ptcp.raw_df[chan][trigger_time+time_delay])       
                    r+=1   


                ptcp.time_delay_signals = {}
                ptcp.time_delay_signals_chans = {}
                for chan in epochs_surrounding:
                    sequences = epochs_surrounding[chan]
                    ptcp.time_delay_signals_chans[chan] = {}
                    for sequence in sequences: # For each 21-length object
                        for time_delay in range(0,len(sequence)):
                            if time_delay in ptcp.time_delay_signals.keys():
                                ptcp.time_delay_signals[time_delay-1].append(sequence[time_delay])
                            else:
                                ptcp.time_delay_signals[time_delay-1] = [sequence[time_delay]]
                            if time_delay in ptcp.time_delay_signals_chans[chan].keys():
                                ptcp.time_delay_signals_chans[chan][time_delay-1].append(sequence[time_delay])
                            else:
                                ptcp.time_delay_signals_chans[chan][time_delay-1] = [sequence[time_delay]]


                ptcp.avgs = []
                for time_delay in ptcp.time_delay_signals.keys():
                    ptcp.avg_ = sum(ptcp.time_delay_signals[time_delay])/len(ptcp.time_delay_signals[time_delay])
                    ptcp.avgs.append(ptcp.avg_)
                    ptcp.std_ = np.std(ptcp.time_delay_signals[time_delay])
                    ptcp.stds = ptcp.avg_/ptcp.std_
                    print(time_delay,ptcp.avg_,ptcp.std_)

                attrs_to_delete = ['raw_df'] # stims
                for attr in attrs_to_delete:
                    try:
                        delattr(ptcp, attrib)
                    except:
                        pass


    def setup_to_run(self):                   
        self.EXP_BASE          = SAVE_DIR # r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Experiments\\'
        os.chdir(self.EXP_BASE)
        self.child_participant_strings, self.adult_participant_strings = set_up_participants() # > set_up_participants found in experiment.py
        print("Number children", len(self.child_participant_strings), "Number adults ", len(self.adult_participant_strings))

