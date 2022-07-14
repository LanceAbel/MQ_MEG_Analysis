from helpers import *
from config import *
from participant_data import *

import sys
import io
# Need this when running in sublime
#sys.stdout = io.TextIOWrapper(sys.stdout.detach(), encoding = 'utf-8')
#sys.stderr = io.TextIOWrapper(sys.stderr.detach(), encoding = 'utf-8')


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
from os.path import exists
import time
import re
import pickle
import hickle
import copy
from pathlib import Path
from IPython.display import clear_output
import inspect


import mne
from mne import io
from mne.stats import permutation_cluster_test
from mne.datasets import sample


class Participant():
    def __init__(self, is_adult_folder, p_id=None, experiment_number=1):
        '''
        is_adult_folder:     desribes which folder to look in, since participants' date are placed in folders specifying whether they are children or not,
        p_id:                is participants' 4-digit code
        experiment_number:   describes which of several .con files to look for, since we take recordings in multiple blocks
        '''
        num_calls = {} # For logger

        self.events_to_tag = events_to_tag
        self.conds_to_compare = conds_to_compare

        self.p_id = p_id
        self.is_adult_folder = is_adult_folder
        self.experiment_number = experiment_number
        self.warnings = [] # Review any issues running the participant after which are not expected to prevent analysis

        self.configure_demographics()                 # Get sex, age etc.
        self.initialise_experiments()                 # Determines which experiments were run on this participant
        
        logger(inspect.stack())  
        
        # Find .con, .mrk, .hsp/.elp files
        self.find_con_file() 
        self.find_mrk_file()
        self.find_headshape_files()
        self.display_file_info()                      # Display which .con, .mrk, .hsp, .elp are being used

        self.get_initials()                           # Finds the participants' initials, and saves them
       
        self.get_raw()                                # Loads the relevant .con file

        logger(inspect.stack())  
            
        
        # Get signals we will be correlating against, annotate raw object
        self.get_surprise()      
        logger(inspect.stack())  
        
        try:
            self.configure_machine()                  # Based on the raw file, works out which machine they were scanned in
        except Exception as e:
            self.warn("Error configuring machine" +str(e))
            self.native_lowpass, self.native_hipass = 200, 0.5
            self.LP_FREQ = min(LP_FREQ,int(self.native_lowpass/1.5))                        
            self.HP_FREQ = max(HP_FREQ,2*self.native_hipass,int(self.native_hipass+0.1))   # 
        
        self.display_raw_info()
        self.applied_z_score = False                  # For appropriately labelling later graphs. This gets set to True if the function is run
        gc.collect()
        


    def report_on_times(self):
        logger(inspect.stack())
        if hasattr(self,'raw'):
            print(" RAW ", self.raw.info)        
        if hasattr(self,'raw_cleaned'):
            print(" RAW CLEANED ", sorted(list(set(self.raw_cleaned.to_data_frame()['time'].values)))[0:100])
            print( "RAW CLEANED", self.raw_cleaned.info['lowpass'], len([x for x in self.raw_cleaned.info['ch_names'] if "MEG" in x]))
        if hasattr(self,'epochs'):
            print( "Epochs ", sorted(list(set(self.epochs.to_data_frame()['time'].values)))[0:100])
            print( "Epochs ", self.epochs.info['lowpass'], len([x for x in self.epochs.info['ch_names'] if "MEG" in x]))    # 10 1ms, 128 ch
        if hasattr(self,'epochs_sound'):
            print( "Epochs Sound", sorted(list(set(self.epochs_sound.to_data_frame()['time'].values)))[0:100])
            print( "Epochs Sound ", self.epochs_sound.info['lowpass'], len([x for x in self.epochs_sound.info['ch_names'] if "MEG" in x]))  # 9 1ms , len -> 0                      
        if hasattr(self,'epochs_ransac'):
            print( "Epochs Ransac", sorted(list(set(self.epochs_ransac.to_data_frame()['time'].values)))[0:100])
            print( "Epochs Ransac", self.epochs_ransac.info['lowpass'], len([x for x in self.epochs_ransac.info['ch_names'] if "MEG" in x]))
        if hasattr(self,'epochs_ransac_autoreject'):
            print( "ERA", sorted(list(set(self.epochs_ransac_autoreject.to_data_frame()['time'].values)))[0:100])
            print( "ERA", self.epochs_ransac_autoreject.info['lowpass'], len([x for x in self.epochs_ransac_autoreject.info['ch_names'] if "MEG" in x]))


    def basic_cleaning(self):
        '''Perform some basic data cleaning steps
        '''
        
        logger(inspect.stack())  
        self.find_bad_channels()    



        logger(inspect.stack())        
        # Cleaning at the channel level   
        self.find_saturations() 
        logger(inspect.stack())        
        #self.find_muscle_artefacts()
        logger(inspect.stack())       

        # Basic filters
        self.bandpass()   

        logger(inspect.stack())   

        self.line_filter()    

        logger(inspect.stack())          

    def initiate_events(self,events_to_tag_dct=events_to_tag, conds_to_compare_dct=conds_to_compare):
        '''Run several steps required for defining events'''
        self.define_events(events_to_tag_dct, conds_to_compare_dct) # Provided above unless otherwise specified

        self.get_surprise()
        logger(inspect.stack())
        self.events_processing()
        logger(inspect.stack())   
          
        
    def define_events(self,events_to_tag_dct,conds_to_compare_dct):
        '''
        Define the events and conditions of interest
        '''
        self.events_to_tag         = events_to_tag_dct
        self.conds_to_compare      = conds_to_compare_dct
        print("Events to tag ", self.events_to_tag )
        print("Conds to compare ", self.conds_to_compare)
        assert sum(self.events_to_tag.keys()) == 1, "Improperly defined events tagging, must sum to 1"
        self.event_to_tag          = self.events_to_tag[1] # Find which event we are tagging
        # e.g. event_to_tag = 'deviants', 

        # The condition choice, then, is one of the tagged events with a specific value (we pick the 2nd one)
        self.cond_A, self.cond_B   = self.conds_to_compare[self.event_to_tag] 
        self.condition_choice      = self.conds_to_compare[self.event_to_tag][1]
        # e.g. if event_to_tag = 'deviants' -> condition_choice = '99',  (cond_A='88',cond_B = '99')
        logger(inspect.stack()  , [self.cond_A, self.cond_B])    
  
        
    def events_processing_define_events(self):
        '''Run several steps required for defining which MEG pulses are events of interest'''
        logger(inspect.stack())    
        if self.event_to_tag == 'frequencies':
            self.define_events_by_tones()
        elif self.event_to_tag == 'deviants':
            self.define_events_by_deviants()        
        elif self.event_to_tag in ['deviants_custom','deviants_specific']:
            self.define_events_by_deviants_custom()  
        elif self.event_to_tag == 'deviants_custom_combine':            
            self.define_events_by_deviants_custom_combine()
        elif 'hgf' in self.event_to_tag or self.event_to_tag=='PS' or self.event_to_tag == 'CS' or self.event_to_tag=='BS':
            print("Will define events by predictor")
            self.define_events_by_predictor()    
               
        logger(inspect.stack())    

    def events_processing_sound(self):
        '''
        i) Find sound events
        ii) Form epochs
        '''
        logger(inspect.stack())  
        self.find_sound_events()  


        logger(inspect.stack())   
        print("Conditions we care about: ", self.conds_we_care_about)   
        print("First events before sound time correction are:", self.events[0:(NUM_TONES+3)])
        self.print_event_codes()
        logger(inspect.stack())    
        #self.plot_events_of_interest()
        #self.find_meg_events()
        # Epochs                 
        self.get_epochs_sound()



        logger(inspect.stack())  
        self.calculate_sound_delays() # Make adjustments to event times
        logger(inspect.stack())  
        self.alter_event_timings()
        print("First events after sound time correction are:", self.events[0:(NUM_TONES+3)])        
        logger(inspect.stack())  
        
        
    def events_processing(self,recalc_events=False):
        '''
        i) Find MEG events
        ii) Define which of the MEG pulses are events of interest
        ii) Associate these with the sound events which shortly follow the MEG pulses
        '''

        logger(inspect.stack())  
        self.find_events(recalc_events=False)

        logger(inspect.stack())  
        self.events_processing_define_events()      
        logger(inspect.stack())  
        #if not hasattr(self, 'sound_delay'): # Don't re-run if done already (e.g. calculating new condition)    # 'channel_num_sound'):  # Has to do so, in order to define epochs_sound
        self.events_processing_sound()
        # else: # Then the pipeline was already run and we know the sound delay, but how that we've re-run the events we must against shift their times
        #     print("Already got sound delay info ")
        #     self.alter_event_timings()


        logger(inspect.stack())  
        self.get_epochs()        # Get new epochs object    


        logger(inspect.stack())  
        self.annotate_epochs()    
        logger(inspect.stack())  


        
    def more_processing(self):
        '''
        Final data processing steps
        '''

        # More data cleaning, will also remove intervals
    
        self.data_cleaning()
        
        logger(inspect.stack())  
        # Calculate signals of interest
        self.signals_of_interest()
        logger(inspect.stack())  
    
        #self.save_participant()

    def data_cleaning(self):
        '''Continue cleaning data with bandpass, line filter, ICA, ransac, autoreject'''

        # Ransac/Autoreject/Prep            
        if RUN_RANSAC:
            self.RANSAC()    
        else: # If we didn't run RANSAC, epochs doesn't change (gets renamed though so it's clear it's the most processed object)
            self.epochs_ransac = self.epochs
        # In case of an error:
        try:
            L = len(self.epochs_ransac)
        except Exception as e:
            self.warn("RANSAC appears to have failed somehow " +str(e) )
            self.epochs_ransac = self.epochs # Handle case where we have skipped RANSAC
        logger(inspect.stack())  
        delattr(self, 'epochs')               
        gc.collect()    


        if RUN_AUTOREJECT:
            self.Autoreject()
        else: # If we didn't run autoreject
            self.epochs_ransac_autoreject = self.epochs_ransac
        # In case of an error:
        try:
            L = len(self.epochs_ransac_autoreject)
        except Exception as e:
            self.warn("Autoreject appears to have failed somehow "+str(e))
            self.epochs_ransac_autoreject = self.epochs_ransac # Handle case where we have skipped RANSAC
        logger(inspect.stack())  
        delattr(self, 'epochs_ransac')           
        gc.collect()             
                        
        # Prep. Any value of building this in?
        if RUN_ICA:
            self.run_ICA()    # Runs after RANSAC and Autoreject 
            
        # Channel alignment can precede data cleaning (no need to do e.g. run Autoreject on discarded channels). No need to run this if not comparing to another group, though.
        # Be aware this will discard data on the adult system, however and that this also changes the standardised signal value on non-discarded channels
        logger(inspect.stack())  

        #self.report_on_times()
        
        #self.align_channels()
        self.align_channels_paul_individual()             
        #self.apply_z_score()                        

    
    def signals_of_interest(self):
        '''Computes signals of interest'''
        
        #logger(inspect.stack())  
        self.epochs_ransac_autoreject.crop(tmin=EPOCH_START_ANALYSIS, tmax=EPOCH_END_ANALYSIS,include_tmax=True) # Can maximise statistical power by restricting our analyses to a known latency range of interest 

        logger(inspect.stack())          
        self.create_evoked_dfs()                   # Average evoked signal over a whole condition
        #self.plot_evoked_dfs()
        logger(inspect.stack())         
        self.create_mmr_dfs()                      # '''Create MMR DFs of a specific condition vs all'''
        logger(inspect.stack())         
        self.create_mmr_difference_dfs(self.cond_A, self.cond_B)        # '''Plot two MMRs (for two conditions) and then the relative MMRs between these two conditions (MMR1 - MMR2)''' 
    
  
    def rerun_with_new_conditions(self,events_to_tag_dct, conds_to_compare_dct):
        '''After defining a new event type, allows you to re-run a different experiment without running all the previous steps'''       
        num_calls = {} # For logger

        logger(inspect.stack())  
        self.initiate_events(events_to_tag_dct,conds_to_compare_dct)
        self.events_to_tag = events_to_tag_dct
        logger(inspect.stack())
        self.downsample()           # We don't need to do this when re-running, since downsampling was done on the first run. But we need to if loading from disk
        self.more_processing()    

        
    def configure_demographics(self):
        '''
        Configure folder and get age/sex/handedness
        '''
        if self.is_adult_folder:
            self.BASE = ADULT_BASE+str(self.p_id)+'\\'
        else:
            self.BASE = CHILD_BASE+str(self.p_id)+'\\'

        try:
            with open(self.BASE+'Age.txt') as f:
                lines = f.readlines()
            self.age = float(lines[0])
        except:
            self.warn("Age file not found for " + str(self.p_id))
            self.age = None       
        try:
            with open(self.BASE+'Sex.txt') as f:
                lines = f.readlines()
            self.sex = str(lines[0])
        except:
            self.warn("Sex file not found for " + str(self.p_id))
            self.sex = None                   
        try:
            with open(self.BASE+'Handedness.txt') as f:
                lines = f.readlines()
            self.handedness = str(lines[0])
        except:
            self.warn("Handedness file not found for " + str(self.p_id))
            self.handedness = None   
            
    def initialise_experiments(self):
        '''Determines, from the named files, which experiment codes were run on this participant
        '''

        self.experiments_ran = {}
        file_names = glob.glob(f"{self.BASE}*.con")
        exp_files = []
        experiments_ran = {}
        for file_name in file_names:
            try:
                match = re.findall(r'_experiment(\d)', file_name)
                if len(list(match)) > 0:
                    self.experiments_ran[match[0]] = None
            except:
                self.warn("Problem with regex on file " +str(file_name))
        print("Experiments ran", self.experiments_ran)

    def find_con_file(self):
        '''Identifies the appropriate .con file which contains the experimental data
        Covers the various corner cases in how previous data was named.
        '''
        
        os.chdir(self.BASE)
        file_names = glob.glob(f"{self.BASE}*.con")  
        self.fname_con = None
        if self.age > 7: # Arbitrary. # Can't look at .con yet!
            try:
                if int(self.p_id) >= CUTOFF_FILE_NAMES: # Then it was data collected by Lance and the experiment is labelled
                    found_experiment = False
                    shortlist_contains = "tspca.con"
                    shortlist = []
                    for file in file_names:
                        if shortlist_contains in file:
                            shortlist.append(file)
                    for f_name in shortlist: # Then there will be one appropriately named file
                        if "experiment"+str(1) in f_name and not found_experiment:
                            self.fname_con = f_name
                            found_experiment = True
                else:
                    try:
                        max_file = max( file_names, key = lambda x: os.stat(x).st_size)         # Largest file (experiments were 15 mins, 10 mins, 10 mins, we want the 15 min one)
                        self.fname_con = max_file   
                    except Exception as e:
                        self.warn("Cannot find .con - 0 " +str(e))
                self.file_identifier = "B1" if "B1" in self.fname_con else None
                self.file_identifier = "B2" if "B2" in self.fname_con else self.file_identifier
                self.file_identifier = "B3" if "B3" in self.fname_con else self.file_identifier
                self.file_identifier = "B4" if "B4" in self.fname_con else self.file_identifier
            except Exception as e:
                self.warn("Cannot find .con - 1 " +str(e))
            if len(self.fname_con) == 0: # Length zero
                self.warn("Cannot find .con - 2 " +str(e))
        else: # ReTHM will likely have been run (it would also likely have found this anyway if the participant's age was old)
            try:
                if int(self.p_id) >= CUTOFF_FILE_NAMES: # Then it was data collected by Lance and the experiment is labelled
                    found_experiment = False
                    shortlist_contains = "rethm_tspca.con"
                    shortlist = []
                    for file in file_names:
                        if shortlist_contains in file:
                            shortlist.append(file)
                    for f_name in shortlist: # Then there will be one appropriately named file
                        if "experiment"+str(1) in f_name and not found_experiment:
                            self.fname_con = f_name
                            found_experiment = True
                else:
                    shortlist_contains = "denoise_rethm.con"
                    shortlist = []
                    for file in file_names:
                        if shortlist_contains in file:
                            shortlist.append(file)
                    if shortlist == []: # Then try a different name
                        shortlist_contains = "denoise.con"
                        for file in file_names:
                            if shortlist_contains in file:
                                shortlist.append(file)
                        if len(shortlist) == 1:
                            self.fname_con = shortlist[0]
                    
                    elif len(shortlist) > 0:
                        if len(shortlist) == 1: # If found it
                            self.fname_con  = shortlist[0]
                        else:
                            self.warn("ALERT: MULTIPLE POTENTIAL .CON FILES ")
                            print("Possible files: ", shortlist)
                            try:
                                max_file = max( shortlist, key = lambda x: os.stat(x).st_size)
                                self.fname_con = max_file

                            except Exception as e:
                                self.warn("Cannot find .con - 3 " +str(e))
                    self.file_identifier = "B1" if "B1" in self.fname_con else None
                    self.file_identifier = "B2" if "B2" in self.fname_con else self.file_identifier
                    self.file_identifier = "B3" if "B3" in self.fname_con else self.file_identifier
                    self.file_identifier = "B4" if "B4" in self.fname_con else self.file_identifier
            except Exception as e:
                self.warn("Cannot find .con - 4 " + str(e))
            if len(self.fname_con) == 0:
                self.warn("Cannot find .con - 5 " + str(e))

        if self.fname_con == None:
            self.warn("Cannot find .con - 6")
       
        #print("CON FILE A: ", self.fname_con, type(self.fname_con))
        self.fname_con = Path(self.fname_con) # Convert string to path
        #print("CON FILE B: ", self.fname_con, type(self.fname_con))

    def warn(self,warning):
        '''Add to warnings'''
        print("%%%%%%%%%%%%%%%% WARNING: "+ str(self.p_id) + " : " + str(warning) + " %%%%%%%%%%%%%%%")
        self.warnings.append(warning)

    def find_mrk_file(self):
        '''
        Finds the marker coil file for this participant in the experiment of interest
        '''
        os.chdir(self.BASE)
        self.fname_mrk = glob.glob(f"{self.BASE}*.mrk")   
        found_mrk = False
        if int(self.p_id) >= CUTOFF_FILE_NAMES:
            for f_name in self.fname_mrk:
                if (str(self.experiment_number)+"_post.mrk" in f_name or str(self.experiment_number)+"_pre.mrk" in f_name) and not found_mrk:
                    self.fname_mrk = f_name   
                    found_mrk = True
        else:
            for f_name in self.fname_mrk:
                if "POST"+self.file_identifier+".mrk" in f_name or "post"+self.file_identifier+".mrk" in f_name or "PRE"+self.file_identifier+".mrk" in f_name or "pre"+self.file_identifier+".mrk" in f_name or "Post"+self.file_identifier+".mrk" in f_name :
                    if not found_mrk:
                        self.fname_mrk = f_name
                        found_mrk = True
                         
        if type(self.fname_mrk) == type(["List!"]):
            for f_name in self.fname_mrk:
                if "Post" in f_name or "Pre" in f_name or "post" in f_name or "POST" in f_name or "pre" in f_name or "PRE" in f_name:
                    if not found_mrk:
                        self.fname_mrk = f_name
                        found_mrk = True
            if not found_mrk:
                self.warn("NO MRK FILE DETECTED WITH EXPECTED NAME. CHOOSING AN ARBITRARY ONE!")
                self.fname_mrk = f_name

    def find_headshape_files(self):
        '''Only one of each file
        '''
        self.fname_elp = glob.glob(f"{self.BASE}*.elp")[0]
        self.fname_hsp = glob.glob(f"{self.BASE}*.hsp")[0]
                
    def display_file_info(self):
        print(self.fname_con)        
        print(self.fname_elp)
        print(self.fname_hsp)
        print(self.fname_mrk)      
        
    def display_raw_info(self): 
        #print(self.raw.info)
        print("Original format :", self.raw.orig_format)
        self.s_freq = self.raw.info['sfreq']
        print("Sampling freq", self.s_freq)
        print("Number of lines ", len(self.raw), self.raw.n_times)
        print("Duration of recording (mins): ", len(self.raw) / self.s_freq / 60 )  

    def get_initials(self):
        self.fname_con = str(self.fname_con)
        self.p_init = re.findall(r'\d\d\d\d_([A-Z][A-Z]).*', self.fname_con)[0]
        print("P_ID: ", self.p_id)
        print("Initials: ", self.p_init)
        with open(self.BASE+str(self.p_id)+"_initials.txt", 'w') as f:
            f.write(self.p_init)
        f.close()
        
    def get_raw(self, trigger_channels_known=False):
        '''Extracts RAW from CON file
        Then works out which system they were scanned in
        '''

        #print(self.fname_con,self.fname_mrk,self.fname_elp,self.fname_hsp,trigger_channels_known,self.is_adult_folder)

        wait_until_memory_free(required_memory = 3.0, max_wait_time_mins = 1) # This requires some memory (set conservatively here as I may run many threads, should only be 500MB)


        is_adult = self.p_id in adult_participant_strings
        if trigger_channels_known:
            print("Trigger channels: ", list(self.trigger_channels))
            stim = self.trigger_channels
        else:
            self.stim_begin = stim_begin_adult if is_adult else stim_begin_child
            stim = list(range(self.stim_begin-1, self.stim_begin-1+NUM_TONES)) #  if not is_adult else list(range(193, 193+NUM_TONES))
        print("Stim is ", stim)
        try:
            self.raw = mne.io.read_raw_kit(input_fname=self.fname_con,  # change depending on file i want
            mrk=self.fname_mrk,
            elp=self.fname_elp,
            hsp=self.fname_hsp,
            stim = stim,
            slope="+",
            stim_code="channel",
            stimthresh=2 if self.is_adult_folder else 1,  # 2 for adults
            preload=True,
            #allow_unknown_format=False,
            verbose=True                     
            )
  
        except Exception as e:
            self.warn("PROBLEM GETTING RAW" + str(e))
            try:
                self.raw = mne.io.read_raw_kit(input_fname=self.fname_con,  # change depending on file i want
                mrk=self.fname_mrk,
                elp=self.fname_elp,
                hsp=self.fname_hsp,
                #stim = list(range(145, 145+NUM_TONES)) if trigger_channels_known == False else self.trigger_channels, #np.arange(193, 200), # Set to child tones (1 less than index for Matlab)
                slope="+",
                stim_code="channel",
                stimthresh=2 if is_adult else 1, # if self.is_adult_folder else 1,  # 2 for adults
                preload=True,
                #allow_unknown_format=False,
                verbose=True
                )
            except Exception as e:
                self.warn("ANOTHER PROBLEM GETTING RAW" + str(e))

        self.fname_con = str(self.fname_con)
        return self.raw

        


    def configure_machine(self):
        '''Configures machine-specific parameters depending on which system they were scanned in
        '''

        if 'self.raw' in locals():
            assert len(self.raw.info.ch_names) in [num_chs_adult, num_chs_child], "Cannot determine which system the person was scanned in"
        else:
            #print("C1: getting RAW")
            self.get_raw(trigger_channels_known=False)
       
        if len(self.raw.info.ch_names) == num_chs_adult:
            self.stim_begin = stim_begin_adult
            self.is_adult_system = True
            self.num_channels = num_meg_chs_adult
            self.audio_channel = audio_channel_adult 
            self.trigger_channels = range(self.stim_begin-1,  self.stim_begin + NUM_TONES-1)  # 193:200 (have to minus 1 from 194:201 as this is Python, not MATLAB)
        elif len(self.raw.info.ch_names) == num_chs_child:
            self.stim_begin = stim_begin_child
            self.is_adult_system = False
            self.num_channels = num_meg_chs_child
            self.audio_channel = audio_channel_child
            self.trigger_channels = range(self.stim_begin-1,  self.stim_begin + NUM_TONES-1) # 145:151 (have to minus 1 from 146:153 as this is Python, not MATLAB)    
        else:
            raise ValueError('Cannot detect the number of channels')
        if self.p_id in SOUND_EXCEPTIONS.keys():
            self.audio_channel = SOUND_EXCEPTIONS[self.p_id]
        
        self.get_raw(trigger_channels_known=True) # Reload given knowledge of trigger channels. Unsure if this is needed
        

        self.native_lowpass = self.raw.info['lowpass'] # Should be 200Hz as set in MEG160 software controlling machine  
        self.native_hipass = self.raw.info['highpass'] # Should be 0.05Hz as set in MEG160 software controlling machine  
        self.LP_FREQ = min(LP_FREQ,int(self.native_lowpass/1.5))                            # min(60,200/1.5 )      -> 60 
        self.HP_FREQ = max(HP_FREQ,2*self.native_hipass,int(self.native_hipass+0.1))        # max(0.5,2*0.5,0.6)    -> 1
        if int(self.p_id) > CUTOFF_FILE_NAMES:
            self.raw.info['experimenter'] = "Lance Abel"
        with open(self.BASE+self.p_id+"_"+self.p_init+"_measurement_time.txt", 'w') as f:
            f.write(str(self.raw.info['meas_date']))
        f.close()
        
    def ch_num_to_name(self,num):
        '''Use zfill instead'''
        if num < 10:
            chan_string = "MEG 00"+str(num)
        elif num < 100:
            chan_string = "MEG 0"+str(num)
        else:
            chan_string = "MEG "+str(num) 
        return chan_string  
    
    def find_bad_channels(self):
        '''Find noisy and flat channels'''
        
        self.raw.info["bads"] = []

        try:
            self.auto_noisy_chs, self.auto_flat_chs, self.auto_scores = find_bad_channels_maxwell(
            self.raw,
            cross_talk=None,
            calibration=None,
            return_scores=True,
            verbose=False,
            ignore_ref=True,
            h_freq=self.LP_FREQ,
            duration=5, # Slice size (interval length)
            )
        except:
            self.warn("Problem finding bad channels..., not using autodetection")
            self.auto_noisy_chs, self.auto_flat_chs, self.auto_scores = [], [], []
            


        print("Noisy channels identified ", self.auto_noisy_chs)  # we should find them!
        print("Flat channels identified", self.auto_flat_chs)   # none for this dataset
        #print("Auto scores", self.auto_scores)       
        
        self.bads = self.raw.info["bads"] + self.auto_noisy_chs + self.auto_flat_chs
        if str(self.p_id) in MANUAL_BAD.keys():
            self.bads = self.bads + MANUAL_BAD[self.p_id]

        bads_list = np.unique(self.bads).tolist()
        self.bads = bads_list
        self.raw.info["bads"] = self.bads
        print("Bad channels including manually labelled ", self.bads)# self.raw.info["bads"]  )


        # Plot bad channels
        logger(inspect.stack(),   self.raw.info.ch_names)
        try:
            self.picks = mne.pick_channels(self.raw.ch_names, include = [], exclude=self.bads) # , exclude=self.bads # regexp='EEG 05.'); shows all channels 50-59
        except Exception as e:
            self.warn("Problem making picks " + str(e))
        logger(inspect.stack(), self.picks)

        #print("Picks (channel nums, with index-1 vs channel name): ", self.picks)
        #picks = mne.pick_channels_regexp(raw.ch_names, regexp='MEG 05.'); # shows all channels 50-59
        #self.raw.plot(order=self.picks, n_channels=10);#len(self.picks));
        #fig = self.raw.plot_psd(tmax=np.inf, fmax=self.LP_FREQ*2, average=True)
        # Plot bad channels within context of good channels
        #self.picks = mne.pick_channels_regexp(self.raw.ch_names, regexp='MEG 05.'); # shows all channels 50-59
        # self.raw.plot(order=self.picks, n_channels=len(self.picks));


        #%% Interpolate bads
        try:
            #self.raw.load_data()
            self.raw = self.raw.interpolate_bads(reset_bads=False,verbose=True) # will clear out raw.info['bads'] 
        except Exception as e:
            self.warn("Problem interpolating bads " + str(e))
        logger(inspect.stack())  
        
        
        
        
    def describe(self):
        '''Just describes existing basic processing results
        '''
        logger(inspect.stack())  
        try:
            print("Bad channels: ", self.bads)
            print("Bad channels: ", self.raw.info['bads'])
        except:
            self.warn("NO BAD CHANS DEFINED")
        try:
            print("Noisy channels identified ", self.auto_noisy_chs)  # we should find them!
        except:
            self.warn("NO NO NOISY CHANS DEFINED")
        try:
            print("Flat channels identified", self.auto_flat_chs)   # none for this dataset
        except:
            self.warn("NO FLAT CHANS DEFINED") 
            
    def find_saturations(self):
        '''Find flat segments (saturations)'''
        print("Bads before saturations ", self.raw.info['bads'])
        annot, bads = mne.preprocessing.annotate_flat(self.raw, bad_percent=3.0, min_duration=0.003, picks=None, verbose=None)
        
        self.raw.info['bads'] = list(set(self.raw.info['bads']+bads))
        self.bads = self.raw.info['bads']
        print("Bads ", bads)
        print("Annot ", annot)
        print("Bads after saturations ", self.raw.info['bads'])

        
    def find_muscle_artefacts(self):
        '''Create annotations for segments that likely contain muscle artifacts'''
        mne.preprocessing.annotate_muscle_zscore(self.raw,ch_type='mag')
        
    def bandpass(self):
        '''Apply bandpass filter'''
        wait_until_memory_free(required_memory = 6, max_wait_time_mins = 1) # This requires some memory (set conservatively here as I may run many threads, should only need a couple GB to run bandpass)

        if hasattr(self,'raw'):        
            self.raw_filt = self.raw.copy().filter(l_freq=self.HP_FREQ, h_freq=self.LP_FREQ)

        #For after running ICA, need to get rid of DC again
        if hasattr(self,'raw_cleaned'):    
            logger(inspect.stack(), "Post-ICA run of bandpass #1")   
            self.raw_cleaned = self.raw_cleaned.copy().filter(l_freq=self.HP_FREQ,h_freq=None)

        if hasattr(self,'epochs_ransac_autoreject'):
            logger(inspect.stack(), "Post-ICA run of bandpass #2")  
            self.epochs_ransac_autoreject = self.epochs_ransac_autoreject.filter(l_freq=self.HP_FREQ,h_freq=None)    

        ## Compare spectra pre and post clean
        # Pre
        # fig = self.raw.plot_psd(tmax=np.inf, fmax=self.LP_FREQ*2, average=True)
        ## Compare spectra pre and post clean
        # Pre
        # fig = self.raw_filt.plot_psd(tmax=np.inf, fmax=self.LP_FREQ*2, average=True)
        
        try_del_attr(self, 'raw')
        gc.collect()      
        

    def line_filter(self):
        wait_until_memory_free(required_memory = 4, max_wait_time_mins = 1) # This requires some memory  (set conservatively here as I may run many threads, should only need 1-2GB to run bandpass)

        '''Line filter option B'''
        self.raw_cleaned = copy.deepcopy(self.raw_filt)
        print("Line filter freqs: ", list(np.arange(LINE_FREQ, LINE_FREQ+LINE_FREQ*8+1, LINE_FREQ)))
        self.raw_cleaned.notch_filter(np.arange(LINE_FREQ, LINE_FREQ+LINE_FREQ*8+1, LINE_FREQ), notch_widths=0.25, n_jobs=NUM_CPUS_LineFilter, filter_length='auto', phase='zero')# , picks=picks) n_jobs=6, # 0.3 keeps the 50Hz line roughly with no bump
        ## Compare spectra pre and post clean
        # Post
        # fig = self.raw_cleaned.plot_psd(tmax=np.inf, fmax=self.LP_FREQ*2, average=True)   
        
        try_del_attr(self, 'raw_filt')
        gc.collect()         

 
    def find_events_from_meg_saved(self):
        '''
        Show the times that tone pules were sent
        '''
        #events = np.genfromtxt(directory+"event_meg.txt", delimiter=',', dtype=None)# skiprows=1)
        #events[:, 0] = self.events[:, 0]
        
   
        file_exists = exists(self.BASE+"event_meg.txt")
        if file_exists:
            f = open(self.BASE+"event_meg.txt","r")
            lines = f.readlines()
            lines_to_add = []
            for line in lines[1:]: # Skip header
                line = line.replace("\n","")
                line = line.split(",")
                newline = [int(line[1]),0,int(line[2]) + self.stim_begin - 2] # Rearrange line and add to the event code (these will be subtracted from) # self.stim_begin -2 = 192 (adults) or 144 (children)  
                newline =np.array(newline)
                lines_to_add.append(newline)
            self.events = np.array(lines_to_add)
        else:
            self.warn("MEG SAVED FILE DOESN'T EXIST")   


    def find_events(self,recalc_events=False):
        '''Renames events from [193, 194, ... ] etc as they appear in the MEG pulses -> integers 1, 2... n'''
        
        logger(inspect.stack())     
#         if recalc_events:
#             raw = self.get_raw()
#             print("Got raw again")
        #else:
        try:
            raw = self.raw_cleaned
            print("Used raw cleaned")
        except Exception as e:
            self.warn("Cannot find raw, re-obtaining " + str(e))
            self.get_raw()
            raw = self.raw
            self.raw_cleaned = self.raw
        
        logger(inspect.stack())  
        self.events = mne.find_events(
            raw,
            output="onset",
            consecutive=False,
            min_duration=0,
            shortest_event=1,# if self.adult else 5,
            mask=None,
            uint_cast=False,
            mask_type="and",
            initial_event=False,
            verbose=None)
        ## Error handling
        if len(self.events) == 0:
            self.warn("Events were empty, looking in saved file")
            self.find_events_from_meg_saved()
            if len(self.events) == 0:
                self.warn("Events were empty, and no saved events file has been found, trying something else....")
                stim = self.trigger_channels
                stim = ["MISC "+str(x-self.num_channels).zfill(3) for x in stim]
                print("Channels we are looking for events on: ", stim)
                self.events = mne.find_events(
                            raw,
                            output="onset",
                            consecutive=False,
                            min_duration=0,
                            shortest_event=1, #if self.adult else 5,
                            mask=None,
                            stim_channel = stim,
                            uint_cast=False,
                            mask_type="and",
                            initial_event=False,
                            verbose=True)                
                if len(self.events) == 0:
                    self.warn("Events are still empty, nothing else I know how to do")

        logger(inspect.stack())   
        self.events_copy = copy.deepcopy(self.events)
        print("First NUM_TONES events: ", self.events_copy[0:NUM_TONES])
        self.events[:, 2] = self.events[:, 2] - (min(self.events[::, 2]) - 1) # Sets them to 1,2,3,4..., NUM_TONES
        logger(inspect.stack())  
        del raw
        gc.collect()        
        
    def define_events_by_tones(self):
        self.events_tones = copy.deepcopy(self.events)
        r = 0
        for event in self.events:
            gap = ALL_TONE_FREQUENCIES[1]-ALL_TONE_FREQUENCIES[0]
            event[-1] = ALL_TONE_FREQUENCIES[0]+gap*(event[-1]-1) # Replace by tone frequency 
            self.events_tones[r] = event[0:]
            r+=1
            
        self.events = self.events_tones        
        self.conds_we_care_about = {str(x): x for x in ALL_TONE_FREQUENCIES}
        
        
    def define_events_by_deviants_common(self):
        '''Reworks events to show the number of repetitions
        Then labels the events as standards or deviants
        '''
        
        self.sequence = []
        self.sequence.append(1)
        j = 2

        for i in range(len(self.events[:, 2]) - 1):
            if self.events[i, 2] == self.events[i + 1, 2]:
                self.sequence.append(j)
                j = j + 1
            else:
                self.sequence.append(1)
                j = 2

        self.repeat_idx = {}
        for j in range(1,MAX_NUM_REPEATS+1):
            self.repeat_idx[j] = np.argwhere(np.array(self.sequence) == j)
            

        self.deviant_idx = np.argwhere(np.array(self.sequence) == 1)
        self.predeviant_idx = (np.argwhere(np.array(self.sequence) == 1) - 1)[1:]  
        print("Pre-deviants: ", self.predeviant_idx[0:5])
        print("Deviants: ", self.deviant_idx[0:5])


        # Define idx_NUM_REPS
        for p in range(1,MAX_NUM_REPEATS+1):
            globals()["self.idx_repeat_"+str(p)] = p
        # Tag the predeviant and deviant along with how many repetitions there had been
        self.idx_predeviant = int(conds_to_compare['deviants'][0])  # At the predeviant
        self.idx_deviant = int(conds_to_compare['deviants'][1])     # Prior to the deviant # Not currently using this to mark anything


        self.new_sequence = np.array(self.sequence)
        for p in range(2,MAX_NUM_REPEATS+1):
            self.new_sequence[self.repeat_idx[j]] = globals()["self.idx_repeat_"+str(p)]
        print("New sequence: ", self.new_sequence[0:(NUM_TONES+3)])   
        
    def define_events_by_deviants(self):
        '''In this world, every event is either a deviant or predeviant.
        Nothing else matters
        '''
        logger(inspect.stack())  
        self.define_events_by_deviants_common()
        logger(inspect.stack())  
        self.new_sequence[self.predeviant_idx] = self.idx_predeviant
        logger(inspect.stack())  
        self.new_sequence[self.deviant_idx] = self.idx_deviant
        logger(inspect.stack())  
        self.events[:, 2] = self.new_sequence  #Override events sequence
        logger(inspect.stack())  
        print("Conditions to compare: ", conds_to_compare)
        self.conds_we_care_about = {conds_to_compare['deviants'][0]: self.idx_predeviant,    conds_to_compare['deviants'][1]: self.idx_deviant}             
        #self.surprise_relevant_series = 'deviants'
        #self.surprise_predictions = self.surprise[self.surprise_relevant_series]
    
    def define_events_by_deviants_custom_generate_deviants(self):
        '''Provide custom tags to events, e.g 996 for a deviant after 6 repetitions
        
        Takes events which were tagged 1,2,3,... based on number of repetitions
        Detects a change from N to 1, labels the N as 88N and the 1 as 99N.
        
        In this world, we have events:
        2,3,4,5,6 # Number of repetitions there have been
        882, 883, 884...88X : Pre-deviants which were the Xth repetition
        991, 992, 993...99X : Deviants which were the Xth repetition
        
        
        '''
        logger(inspect.stack())  
        self.define_events_by_deviants_common()

        logger(inspect.stack(), self.events)
        self.events[:, 2] = self.new_sequence #Override events sequence with 88/99 for pre-deviants and deviants
        
        # Type 1 - marks how many repetitions before the deviant
        last_event = 1
        r = 0
        self.all_event_codes = []
        for event in self.events:
            self.all_event_codes.append(event[2])

        # Custom deviant tagging
        r = 0
        for r in range(0,len(self.all_event_codes)):           
            if r>=1:
                if self.all_event_codes[r] > self.events[r-1][2]:   # A standard
                    last_event = self.all_event_codes[r]
                else: # Ordinary deviant, or double deviant
                    if self.all_event_codes[r] < self.events[r-1][2]:  # Ordinary deviant
                        self.all_event_codes[r] = int(conds_to_compare['deviants'][1] +str(last_event))
                        self.all_event_codes[r-1] = int(conds_to_compare['deviants'][0]+str(last_event))               
                        last_event = 1    
                    elif self.all_event_codes[r] == 1: # 2nd deviant in a row
                        self.all_event_codes[r] = int(conds_to_compare['deviants'][1]+"1") # 991
                        last_event = 1    

    
    def define_events_by_deviants_custom_generate_events(self):    
        ''' Replace event codes with the custom deviants codes
        These track how many repetitions there were before a deviant (i.e. 886=6th repetition, 996=deviant after 6 repetitions)
        ''' 

        p = 0
        for event in self.events:
            self.events[p][2] = self.all_event_codes[p]
            p+=1
        self.conds_we_care_about = {}
        # We still define repetitions:  
        for j in range(2,MAX_NUM_REPEATS): # 1 is a deviant so is always labelled as such. # Event "7" becomes event 887 always, as it is always a pre-deviant
            self.conds_we_care_about[str(j)] = globals()["self.idx_repeat_"+str(j)]
        # Predeviant
        for p in range(882,888):
            self.conds_we_care_about[str(p)] = p 
        # Deviant
        for p in range(991,998):
            self.conds_we_care_about[str(p)] = p  
    
    def define_events_by_deviants_custom(self):
        '''
        Define events as deviant or predeviant with custom tags based on how many repetitions there have been before the deviant
        '''
        
        logger(inspect.stack())  
        self.define_events_by_deviants_custom_generate_deviants()
        logger(inspect.stack())  
        self.define_events_by_deviants_custom_generate_events()
        logger(inspect.stack())  

    def define_events_by_deviants_custom_combine(self):
        '''Combines several categories of custom deviants
        
        These categories now represent:
        category B: deviants after 5-7 repeats ('999') vs
        category A: deviants after 1-4 repeats ('991)
        '''
        logger(inspect.stack())  
        self.define_events_by_deviants_custom_generate_deviants()
        logger(inspect.stack())  
        self.define_events_by_deviants_custom_generate_events() #added
        print("First 10 events ", self.events[0:10])
        logger(inspect.stack())  
        # Make replacements
        p = 0
        for event in self.events:
            if self.all_event_codes[p] >= 995:
                self.events[p][2] = 999
            elif self.all_event_codes[p] >= 991:
                self.events[p][2] = 991
            p+=1
        
        self.conds_we_care_about = {}
        # Deviants after few repetitions
        self.conds_we_care_about[str(991)] = 991
        # Deviants after many repetitions
        self.conds_we_care_about[str(999)] = 999           

   
    def define_events_by_predictor(self):
        '''
        Load the predictor (e.g. HGF)
        In this HGF, every event's surprise is defined by an independent predictor
        '''

        self.define_events_by_deviants_common()

        ## Below is verbose, fix.
        ## Layer 2
        if self.events_to_tag[1] in ["hgf_pe2"]: # ["hgf_pe2", "hgf_pe2_extremes"]:
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PE2"             # Vanilla HGF PE2 (hgf_whatworld)
        elif self.events_to_tag[1] == "hgf_pe2_mod":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PE2_mod"         # HGF Modified PE2 - The importance of tone frequency is hacked onto the end to modify the PE2 and PWPE2            
        elif self.events_to_tag[1] == "hgf_pe2_mod_baked":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PE2_mod_baked"   # HGF Integrated PE2 - The importance of tone frequency is based into the modified HGF model            
        
        elif self.events_to_tag[1] == "hgf_pwpe2":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PWPE2"           # Vanilla HGF PWPE2 (hgf_whatworld)
        elif self.events_to_tag[1] == "hgf_pwpe2_mod":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PWPE2_mod"       # HGF Modified PWPE2    
        elif self.events_to_tag[1] == "hgf_pwpe2_mod_baked":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PWPE2_mod_baked" # HGF Integrated PWPE2     
            
        ## Layer 3
        elif self.events_to_tag[1] == "hgf_pe3":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PE3"             # Vanilla HGF PE3 (hgf_whatworld)
        elif self.events_to_tag[1] == "hgf_pe3_mod":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PE3_mod"         # HGF Modified PE3               
        elif self.events_to_tag[1] == "hgf_pe3_mod_baked":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PE3_mod_baked"   # HGF Integrated PE3 
            
        elif self.events_to_tag[1] == "hgf_pwpe3":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PWPE3"           # Vanilla HGF PWPE3 (hgf_whatworld)
        elif self.events_to_tag[1] == "hgf_pwpe3_mod":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PWPE3_mod"       # HGF Modified PWPE3   
        elif self.events_to_tag[1] == "hgf_pwpe3_mod_baked":
            self.surprise_relevant_series = "HGF_"+str(RELEVANT_LENGTH_PREDICTOR) + "_PWPE3_mod_baked" # HGF Integrated PWPE3  
        ## Different models of surprise
        elif self.events_to_tag[1] in ["BS","PS","CS"]:
            self.surprise_relevant_series = self.events_to_tag[1]
        print("Relevant predictor series ", self.surprise_relevant_series)
            
        ## Get the surprise values for the specific predictor, and the cutoff value that gives a "high" level of surprise. These should be about as frequent as a deviant relative to a pre-deviant
        self.surprise_predictions = self.surprise[self.surprise_relevant_series]
        # print("Surprise predictions " , self.surprise_predictions)
        ## Get the cutoff in value for surprise predictions which partitions it (e.g. into low vs high surprise)
        self.surprise_cutoff = float(self.surprise_cutoffs[self.surprise_relevant_series])
        if "hgf" in self.events_to_tag[1] or "CS" in self.events_to_tag[1] or "PS" in self.events_to_tag[1] or "BS" in self.events_to_tag[1]:
            if CONT_PREDICTOR_USE_DIFF_CUTOFF: # Then use a different low and high cutoff
                self.surprise_cutoff_low = float(self.surprise_cutoffs_low[self.surprise_relevant_series])
            else:
                self.surprise_cutoff_low = self.surprise_cutoff
        else:
            self.surprise_cutoff_low = self.surprise_cutoff
        print("Surprise cutoff ", self.surprise_cutoff)            
        print("Surprise cutoff low ", self.surprise_cutoff_low)
        
        ## Label the event with reference to the predicted surprise value
        self.surprise_category = []
        for predicted_surprise in self.surprise_predictions:
            if predicted_surprise >= self.surprise_cutoff:
                self.surprise_category.append(int(conds_to_compare[self.events_to_tag[1]][1])) # Tag for high surprise
            elif predicted_surprise < self.surprise_cutoff_low:
                self.surprise_category.append(int(conds_to_compare[self.events_to_tag[1]][0])) # Tag for low surprise              
            else: # (i.e. self.surprise_cutoff != self.surprise_cutoff_low)
                self.surprise_category.append(666)          # No man's land. Neither high nor low preedicted surprise  # Ignore these later, if they exist
        ## Apply the cutoff by labelling the events according to the predictor
        print("Length of events: ", len(self.events), "Events: ", self.events[0:10])
        print("Num in category: ", len(self.surprise_category))
        print("Category classifications: ", set(self.surprise_category))
        min_length = min(len(self.events),len(self.surprise_category))
        self.events = self.events[0:min_length]
        self.events[:, 2] = self.surprise_category[0:min_length]  #Override events sequence (need at least as many predictions as there are tones!)

        self.conds_we_care_about = {}
        for event_label in conds_to_compare[self.events_to_tag[1]]:
            self.conds_we_care_about[event_label] = int(event_label)
         
    def print_event_codes(self):
        '''Print a list of all event codes'''
        all_evts = []
        for evt in self.events:
            all_evts.append(evt[2])
        all_evts = set(all_evts)
        print("All event codes: ", list(all_evts))
        
 
    def plot_events_of_interest(self):
        # try:
        #     fig = mne.viz.plot_events(self.events, event_id=self.conds_we_care_about, sfreq=self.raw_cleaned.info['sfreq'],
        #                               first_samp=self.begin_time*self.s_freq)
        # except Exception as e:
        #     print("Error probably not found: did you sample a large enough number of events")
        #     print("Exception", str(e))

        #%% Find sound onset delay relative to portcode output (detected event time)
        picks = mne.pick_channels(self.raw_cleaned.info["ch_names"], include=[self.channel_name_sound])  
        

    def find_meg_events(self):
        '''
        Show the times that tone pules were sent
        '''
        
        f = open(self.BASE+"event_meg.txt","r")  # Need to run MATLAB script first which generates this
                                            # Path: E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Roving_MMN_2020-master\Roving_MMN_2020-master\preprocessing_scripts_ME125_phase2_YSun\
                                            # Filename: ME125_Phase2_Step5.m
        lines = f.readlines()
        self.tone_sequence = []
        for line in lines[1:]: # Skip header
            line = line.replace("\n","")
            tone_num = line[-1]
            self.tone_sequence.append(tone_num)
        print("Tone sequence :", self.tone_sequence[0:MAX_NUM_REPEATS])            
            
    def find_sound_events(self):
        '''
        Plots the sound events
        '''
        logger(inspect.stack())  
        #self.channel_name_sound = 'MISC 007'
        self.channel_num_sound = self.audio_channel - self.num_channels   # 167-160 = 7 = "MISC 007", or 135-125 = "MISC 010"       There are 3x reference sensors 126,127,128 in kids (161,162,163 in adults), then a few others
        self.channel_name_sound = "MISC "+str(self.channel_num_sound).zfill(3)
        logger(inspect.stack())  
        # Check sound coming through
        self.begin_time    = 45 # secs into recording. Sometimes, takes a while to start the script
        self.time_to_plot  = 30  # Enough to see the sound plays on *every* trigger
        self.end_time = self.begin_time + self.time_to_plot
        rows_to_plot = range(int(self.begin_time*self.s_freq), int(self.end_time*self.s_freq))
        logger(inspect.stack())  
        try:
            self.create_df()
            self.times = self.df.index.values
            # plt.plot(self.df[self.channel_name_sound][self.times[rows_to_plot]])
            # plt.legend(["Sound channel data: "], bbox_to_anchor=(0.75, 1.15), ncol=2)
            # plt.show()
        except Exception as e:
            self.warn("Could not extract times from dataframe " + str(e))
            self.times = [event[0] for event in self.events]        
        logger(inspect.stack())  
        try_del_attr(self,'df')
        gc.collect()
    

    def get_epochs_sound(self):
        '''We only care about the sounds when they are for conditions that we care about (that could be for every tone, though)
        This ignores all channels other than the one containing sound
        '''
        self.epochs_sound = mne.Epochs(
                                        self.raw_cleaned,
                                        self.events,
                                        event_id=self.conds_we_care_about, # Which event tags to count as events
                                        tmin=EPOCH_START,
                                        tmax=EPOCH_END,
                                        #decim=5, # subsampling ; 1000/decim = new frequency
                                        #picks=self.picks,
                                        baseline=(None,None),
                                        preload=True,
                                        )
        
    def get_sound_times(self):
        '''Load sound times from file'''
        try:
            with open(self.BASE+'event_sounds.txt') as f:
                lines = f.readlines()
            self.sound_delay = float(lines[0])*1/SEC_TO_MS
            print("New sound delay %s ms "%(str(self.sound_delay*SEC_TO_MS)))
        except:
            self.warn("Cannot find sound events file")

    def get_saved_sound_delay(self):
        '''Gets median sound delay from saved file if this isn't calculated successfully
        This analysis takes place in compare_timing.py
        '''
        try:
            with open(self.BASE+'Sound Delay.txt') as f:
                lines = f.readlines()
            self.sound_delay = float(lines[0])*1/SEC_TO_MS
            print("New sound delay %s ms "%(str(self.sound_delay*SEC_TO_MS)))
        except:
            self.warn("Cannot find saved sound delay file")

    def calculate_sound_delays_from_file(self):
        '''
        Populate a list of sound delays used for correcting the event time
        '''
        try:
            # Load sounds
            with open(self.BASE+'event_sounds.txt') as f:
                events_sound = f.readlines()
            events_sound = [event_sound.replace("\n","") for event_sound in events_sound]
            times_sound = events_sound[1:]
            times_sound = [time_sound.split(",") for time_sound in times_sound]
            
            # Load MEG pulses
            with open(self.BASE+'event_meg.txt') as f:
                events_meg = f.readlines()
            events_meg = [event_meg.replace("\n","") for event_meg in events_meg]        
            times_meg = events_meg[1:]        
            times_meg = [time_meg.split(",") for time_meg in times_meg]
            
            assert len(times_meg)==len(times_sound), "There are different numbers of events and sounds, cannot replace event times using sound file"
            diffs = [float(times_meg[r][1])-float(times_sound[r][1]) for r in range(0,len(times_meg))]
            median_diff = np.median(diffs)
            self.sound_delay = -median_diff*1/SEC_TO_MS # convert ms to secs
            print("Median diff", median_diff, "std", np.std(diffs)) 
            
            # Replace aberrant differences with the median
            error_robust_diffs = []
            error_samples = 0
            for diff in diffs:
                if abs(diff-median_diff) < 10:
                    error_robust_diffs.append(diff)
                else:
                    error_robust_diffs.append(median_diff)
                    error_samples+=1
            if error_samples > 0:
                self.warn("Replaced %s out of %s sound events with the median delay time, the rest were replaced based on actual delay time.. "%(str(error_samples),str(len(diffs))))
            
            r = 0
            for evt in self.events:
                self.events[r][0] = round(evt[0]-error_robust_diffs[r],0) #ms
                r+=1
            print("First 5 events : ", self.events[0:5])    
            self.times = [event[0] for event in self.events]
            logger(inspect.stack(), "Successfully replaced all event times based on actual sound delays") 
            self.sound_delay_from_file = True
        except Exception as e:
            self.warn("Could not use sound events file" + str(e))


    def calculate_sound_delays(self):
        '''
        Adjust the timestamps of the events for the sound delay
        Plot the typical timecourse of the sound during the epoch
        '''

        self.calculate_sound_delays_from_file()
        ## If we couldn't work out sound delay from sound events file, then try to use the sound wave itself to do so...
        if not hasattr(self,'sound_delay_from_file'):
            try:            
                sound_av = np.abs(self.epochs_sound.pick_channels([self.channel_name_sound]).get_data())  # account for the fact that undersampled sound varies in polarity
                #sound_av = self.df[self.channel_name_sound].values
                sound_av = np.mean(np.squeeze(sound_av, axis=None), axis=0)
                sound_av = sound_av - np.mean( sound_av[0 : int(np.argwhere(self.epochs_sound.times == 0))])  # Remove any offset

                max_prestim = np.max(sound_av[0 : int(np.argwhere(self.epochs_sound.times == 0))])
                mean_prestim = np.mean(sound_av[0 : int(np.argwhere(self.epochs_sound.times == 0))])
                self.sound_delay = np.min(self.epochs_sound.times[sound_av > 3 * max_prestim])
            except Exception as e:
                self.warn("No Usable Audio Pulses found: " +str(e) + " using saved sound delay ")
                if int(self.p_id) >= SOUND_DELAY_2_BEGINS:
                    self.sound_delay = SOUND_DELAY_2
                else:
                    self.sound_delay = SOUND_DELAY_1
                self.get_saved_sound_delay()
            # plt.plot(self.epochs_sound.times, sound_av[:], self.sound_delay, mean_prestim, "g.")
            # plt.legend(["Epochs sound "], bbox_to_anchor=(0.75, 1.15), ncol=2)
            # plt.show()      
         
            del sound_av
            gc.collect()
        
    def alter_event_timings(self):
        '''Adjust events for sound onset delay '''
        self.events[:, 0] = self.events[:, 0] + int(np.round(self.sound_delay * SEC_TO_MS))    
        self.times = [event[0] for event in self.events]        
        print("Sound delay %s ms "%(str(self.sound_delay*SEC_TO_MS)))
    
    def get_epochs(self,downsample=False):
   
        '''Epoching - wide for checking that timings are correct'''
        self.epochs = mne.Epochs(
            self.raw_cleaned,
            self.events,
            event_id=self.conds_we_care_about,
            tmin=EPOCH_START,
            tmax=EPOCH_END,
            preload=True,
            baseline=(None,None))

        # # Resample if necessary
        # if downsample==True:
        #     print("DOWNSAMPLING...")
        #     self.epochs = self.epochs.resample(DOWNSAMPLE_FREQ)

        #%% Choose everything - unecessary here but can be used to filter classes
        self.picks = mne.pick_types(self.epochs.info, meg="mag", eeg=False, stim=False, eog=False, include=[], exclude=[])

    def annotate_epochs(self):
        '''Adds labelling to all epochs
        '''
        all_annotations = []
        sound_annot = mne.Annotations(onset=self.events[:, 0]/SEC_TO_MS,  # in seconds
                                       duration=[TONE_DURATION] * len(self.events[:, 0]),  # in seconds, too
                                       description=['sound'] * len(self.events[:, 0]))
        self.raw_cleaned.set_annotations(sound_annot)
        #self.epochs.set_annotations(sound_annot)
        
    def plot_ecg(self):
        '''Plot channels identified as containing heartbeats'''
        ecg_chans_to_plot = list(range(51,58)) + list(range(96,102))+list(range(112,119)) # Just as a sample
        self.picks = mne.pick_channels(self.raw_cleaned.ch_names, include=[], exclude=self.bads); # regexp='EEG 05.'); shows all channels 50-59
        self.raw_cleaned.plot(order=ecg_chans_to_plot, n_channels=20);#len(self.picks));
        
    def plot_eog(self):
        ''' Plot channels identified as containing blinks'''
        eog_chans_to_plot = list(range(36,45)) + list(range(46,48))+list(range(101,112))
        self.picks = mne.pick_channels(self.raw_cleaned.ch_names, include=[], exclude=self.bads); # regexp='EEG 05.'); shows all channels 50-59
        self.raw_cleaned.plot(order=eog_chans_to_plot, n_channels=20);#len(self.picks));
    
    def save_fif(self):
        print(self.raw[-1, :]) # last column is stim channel
        des = "_lp_hp_badsremoved" # Description of operations done thus far
        self.raw.save(fname+des+'.fif', picks=None, buffer_size_sec=None, drop_small_buffer=False, proj=False,fmt='single',overwrite=True, split_size='2GB', split_naming='neuromag', verbose=None)
        
    def create_df(self):
        try:
            self.df = self.raw_cleaned.to_data_frame()
        except Exception as e:
            self.warn("Exception making dataframe " + str(e))

    def downsample(self):
        '''Resample if required'''
        # Note this is to samples/sec
        down_sFreq = 200 # int(self.LP_FREQ*0.5) + 10         # int(60*0.5) + 10 -> 40
        self.raw_cleaned.resample(down_sFreq)           
        self.epochs.resample(down_sFreq) # Speeds up next steps, such as RANSAC/Autoreject
        self.s_freq = self.raw_cleaned.info['sfreq']

    def RANSAC(self):
        '''
        Apply Ransac cleaning (Jas, M., Engemann, D. A., Bekhti et al, 2017)
        '''

        self.ransac = Ransac(verbose='progressbar', picks=self.picks, n_jobs=NUM_CPUS_Ransac)  #  not sure if this (parallel jobs) actually results in a speedup?  # verbose = True, 
        self.epochs_ransac = self.ransac.fit_transform(self.epochs) 
        self.reject_threshold = get_rejection_threshold(self.epochs_ransac)
        # Reject threshold
        print("The rejection dictionary is %s" % self.reject_threshold)        
        self.epochs_ransac.drop_bad(reject=self.reject_threshold,verbose=False)

        self.print_conditions()
        
        #print("Plot RANSAC drop log")
        #self.epochs_ransac.plot_drop_log()#verbose=False);

        ## These events can be equalised, but they would already have very similar numbers of events...
        # if self.events_to_tag[1] in ["deviants_custom_combine", "deviants_specific"]:
        #     # Equalise events
        #     print("Equalising event counts")
        #     self.epochs_ransac.equalize_event_counts(self.conds_we_care_about) 

        # Plot a few channels
        # for channel_num in range(1,5):
        #     self.epochs_ransac.plot_image(picks=['MEG 00'+str(channel_num)]) # Both methods have a picks parameter for subselecting which channel(s) to return; raw.get_data() has additional parameters for restricting the time domain
        #self.time_freq_morlet()
    
    
    def time_freq_morlet(self):
        pass
        #frequencies = np.arange(7, 30, 3)
        
#         power = mne.time_frequency.tfr_morlet(self.epochs_ransac, n_cycles=2, return_itc=False,
#                                               freqs=frequencies, decim=3)
#         power.plot(['MEG 001'])    
 
    def Autoreject(self):
        
        '''Run AUTOREJECT (Jas, M., Engemann, D. A., Bekhti et al, 2017)'''

        time_ = time.time()
        consensus_percs = np.linspace(0, 1.0, 11)
        # The values to try for percentage of channels that must agree as a fraction of the total number of channels. This sets kappa/Q. If None, defaults to np.linspace(0, 1.0, 11)
        # These are the kappa (κ) values that autoreject will try (see the autoreject paper) for more information on these parameters).
        # e.g. if 0.4, if >40% of channels exceed the channel threshold then the epoch is dropped

        n_interpolates = np.array([1, 4, 32])
        # The values to try for the number of channels for which to interpolate. This is rho. If None, defaults to np.array([1, 4, 32])
        # These are the rho (ρ) values that we would like autoreject to try 
        # e.g. if 7, then, for channels that are not discarded, the worst n_interpolates channels are interpolated

        
        # We need to temporarily clear the "bads" list, as autoreject doesn't run on these (and we want it to as much as the channels from which it was interpolated)
        # Start added
        bads_temp = copy.deepcopy(self.epochs_ransac.info["bads"][0:])
        bads_temp2 = copy.deepcopy(self.bads[0:])
        picks_temp = self.picks[0:]
        self.epochs_ransac.info["bads"] = []
        self.bads = []
        #self.picks = mne.pick_channels(self.epochs_ransac.ch_names, include=[], exclude=self.bads)
        # End added
        self.ar = AutoReject(n_interpolates, consensus_percs, picks=self.picks, thresh_method='random_search', random_state=42, verbose=False, n_jobs=NUM_CPUS_Autoreject);
        start_time = time.time()
        
        
        # Detect which epochs object to use
        # if not self.epochs_ransac in locals() and not self.epochs_ransac in globals(): 
        try:
            L = len(self.epochs_ransac)
        except:
            self.epochs_ransac = self.epochs_sound # Handle case where we have skipped RANSAC
            self.warn("Ransac epochs not found, using uncleaned sound epochs")  
        try_del_attr(self, 'epochs_sound')
        gc.collect()          

        # Run autoreject
        try:
            self.ar.fit(self.epochs_ransac);
        except Exception as e:
            self.warn("Error running AUTOREJECT " +str(e))   
        print("Running autoreject1 took: ", time.time()-start_time, " seconds")   

        
        #print("Plotting before autoreject")
        #self.epochs_ransac.plot(n_epochs=10);
        self.epochs_ransac_autoreject, self.reject_log = self.ar.transform(self.epochs_ransac, return_log=True)
        print("Running autoreject2 took: ", time.time()-start_time, " seconds") 
        

        ## Plot after autoreject
        #print("Plotting after autoreject")
        #self.epochs_ransac_autoreject.plot(n_epochs=10);       


        ## Visualize the dropped epochs (must refer to old object)
        #print("METHOD 1: Visualising dropped epochs")
        # try:
        #     self.epochs_ransac[self.reject_log.bad_epochs].plot(scalings=dict(mag=4000e-15))
        # except:
            # self.warn("Error visualising dropped epochs")   

        ## Plot auto-rejections
        try:
            pass
            #print("METHOD 2: Plotting reject log")
            # self.reject_log.plot()    
            # Good data (light red)
            # Bad segment but not to be interpolated (medium dark red)
            # Bad segment to be interpolated (dark red)
        except:
            self.warn("Error plotting auto-reject log")
        
        try: 
            pass
            #print("METHOD 3: Plotting epochs auto-rejected")
            # self.reject_log.plot_epochs(self.epochs_ransac);    
        except:
            self.warn("Error plotting epochs auto-rejected")
            
        ## Display autoreject actions
        self.reject_threshold = get_rejection_threshold(self.epochs_ransac) # Outlier?
        print("Auto-reject thresholds : ", self.reject_threshold)

        ## Get bad epochs
        self.reject_indices = [i for i in range(0,len(self.reject_log.bad_epochs)) if self.reject_log.bad_epochs[i] == True]
        print("Auto-reject indices :", self.reject_indices)  
        
        ## Restore bads
        self.epochs_ransac.info["bads"] = bads_temp[0:]
        self.bads = bads_temp2[0:]
        #self.picks = picks_temp
        
    def run_ICA(self):
        '''Cleaning data using ICA; removes components unrelated to brain activity (cardiac, ocular)'''

        logger(inspect.stack())  

        self.ica = ICA(n_components=ICA_NUM_COMPONENTS, max_iter="auto", method="fastica", random_state=97,verbose=False)
        self.ica.exclude = []
        self.ica.excludes = {}  # For record-keeping purposes
        logger(inspect.stack())  

        reject = dict(mag=5e-12)    #, grad=4000e-13) # Avoid fitting ICA on crazy environmental artifacts that would dominate the variance and decomposition
        #reject = dict(mag=3e-12)   # Paul suggestion was 5e-13, as it's not really a magnetomer but a gradiometer. This produced error RuntimeError('No clean segment found. Please consider updating your rejection thresholds.',)]]
                                    # 5e-12 was the smallest I could get it without producing an error
        
        picks_meg = mne.pick_types(self.raw_cleaned.info, meg=True, eeg=False, eog=False,stim=False, exclude='bads')


        ## Could fit to entire raw data or to epochs, if not all of the data belonged to an epoch, this would then get less data
        self.ica.fit(self.raw_cleaned,picks=picks_meg, reject=reject)
        #self.ica.fit(self.epochs_ransac_autoreject,picks=picks_meg,reject=reject)
        self.ica_copy = copy.deepcopy(self.ica)

        ## Plots
        #self.ica.plot_components();
        logger(inspect.stack())  
        #self.ica.plot_properties(self.raw_cleaned)
        #self.ica.plot_properties(self.epochs_ransac_autoreject,psd_args={'fmax': 35.})


        self.ica_manual_exclude_list = [] # eval(input("List of ICA components to manually remove: ")) # Input is string
        ## All ICA
        #if self.ica_manual_exclude_list in locals().keys() or self.ica_manual_exclude_list in globals().keys():
        self.ica.exclude = self.ica.exclude + self.ica_manual_exclude_list
        self.ica.exclude = list(set(self.ica.exclude))
        logger(inspect.stack())  


        #print("Plotting ICA before removal")
        #self.ica.plot_sources(self.epochs_ransac_autoreject, show_scrollbars=False);
        #self.ica.plot_sources(self.raw_cleaned, show_scrollbars=False);    
        logger(inspect.stack())  
        self.ICA_ECG()
        logger(inspect.stack())  
        self.ICA_EOG()
        logger(inspect.stack())  
        self.apply_ICA()
        self.bandpass() # Run again as ICA can introduce DC shifts
        logger(inspect.stack())  
        #self.ica.plot_sources(self.raw_cleaned, show_scrollbars=False);  
        if not SAVE_ICA:
            del self.ica
        gc.collect()


    def ICA_EOG(self, min_channels_ica=MIN_CHANNELS_ICA):
        '''Identify EOG components'''
        print("ICA EOG")
        self.all_eog_scores = {}
        self.all_eog_indices = {}
        logger(inspect.stack())  

        for p in range(1,self.num_channels): # Typical ECG channels if head is centred in child system: self.list(range(51,58)) + list(range(96,102))+list(range(112,119)):
            channel_name = "MEG "+str(p).zfill(3)
            self.eog_indices, self.eog_scores = self.ica.find_bads_eog(self.raw_cleaned, ch_name=channel_name, measure='correlation', threshold = ICA_THRESHOLD_EOG ) #   # find which ICs match the EOG pattern 
            # Probably should run on continuous data as EOG events aren't that frequent/may not appear in most epochs
            
            #print("Channel %s , EOG indices: %s "%(str(channel_name),self.eog_indices))
            #self.eog_epochs.plot_image(combine='mean')        
            self.all_eog_scores[p] = self.eog_scores
            self.all_eog_indices[p] = self.eog_indices
            gc.collect()          

        logger(inspect.stack())  
        try:
            self.eog_bads = []
            self.eog_bad_ct = {}
            for index_all_channels in self.all_eog_indices.values():
                for ica_index in index_all_channels:
                    self.eog_bads.append(ica_index)
                    if ica_index in self.eog_bad_ct.keys():
                        self.eog_bad_ct[ica_index] = self.eog_bad_ct[ica_index] +1
                    else:
                        self.eog_bad_ct[ica_index] = 1
            self.eog_bads = list(set(self.eog_bads))
            print("ALL EOG bads ", self.eog_bads)
            print("All EOG bad Counts: ", self.eog_bad_ct)

            # Get two most frequently identified ICA indices to exclude
            max_key = max(self.eog_bad_ct, key= self.eog_bad_ct.get) if len( self.eog_bad_ct) > 0 else None
            second_largest_key = None
            ica_counts = self.eog_bad_ct.values()
            if len(self.eog_bad_ct) > 1:
                second_largest = sorted(ica_counts)[max(-2,-len( self.eog_bad_ct))]
            else:
                second_largest = None

            for count_num in ica_counts:
                for ica_index in self.eog_bad_ct.keys():
                    if self.eog_bad_ct[ica_index] == second_largest:
                        second_largest_key = ica_index
            print("Max ", max_key, " 2nd: ", second_largest_key)
            self.ica.excludes['EOG'] = [max_key,second_largest_key]
            # Exclude the two most frequent ICA components, if there are enough channels identifying this component
            if second_largest_key != None and self.eog_bad_ct[second_largest_key] >= min_channels_ica:
                self.ica.exclude = self.ica.exclude + [second_largest_key]
            if max_key!=None and self.eog_bad_ct[max_key] >= min_channels_ica:
                self.ica.exclude = self.ica.exclude + [max_key]
        except:
            print("Issue getting eog indices")


        # Specifically for EOG/ECG
        if self.eog_bads:
            logger(inspect.stack())  
            # try:
            #     self.ica.plot_sources(self.raw_cleaned, picks=self.eog_indices)
            # except:
            #     print("Problem plotting ICA sources, were there EOG indices?")
            # try:
            #     self.ica.plot_properties(self.raw_cleaned, picks=self.eog_indices, psd_args={'fmax': 35.})        
            # except:
            #     print("Problem plotting ICA properties, were there EOG indices?")
            #self.ica.plot_sources(self.epochs_ransac_autoreject, picks=self.eog_indices);
            #self.ica.plot_properties(self.epochs_ransac_autoreject, picks=self.eog_indices, psd_args={'fmax': 35.});
            try:
                self.eog_evoked = mne.preprocessing.create_eog_epochs(self.raw_cleaned).average()  ## Arbitrary channel (should not be like this). Why do I need to add one but not for ECG?
                #self.eog_evoked.plot_image(combine='mean')
                #self.eog_evoked.plot(gfp='only', spatial_colors=True, ylim=dict(eeg=[-12, 12]))            
            except:
                self.warn("Issue getting EOG evoked")  
            # barplot of ICA component "EOG match" scores
            #self.ica.plot_scores(self.eog_scores);

    def ICA_ECG(self, min_channels_ica=MIN_CHANNELS_ICA):
        '''Identify ECG components'''
        print("ICA ECG")

        self.all_ecg_indices = {}
        self.all_ecg_scores  = {}
        self.all_ecg_events  = {}
        self.all_ecg_epochs  = {}
        logger(inspect.stack())  
        for p in range(1,self.num_channels): # Typical ECG: self.list(range(51,58)) + list(range(96,102))+list(range(112,119)): #range(1,10):#25): 
            channel_name = "MEG "+str(p).zfill(3)
            self.ecg_indices, self.ecg_scores = self.ica.find_bads_ecg(self.raw_cleaned, ch_name=channel_name, method='correlation', threshold=ICA_THRESHOLD_ECG ) # method="ctps"  # find which ICs match the ECG pattern     
            #self.ecg_indices, self.ecg_scores = self.ica.find_bads_ecg(self.raw_cleaned, method="ctps" ) #   # find which ICs match the ECG pattern             
            print("Channel %s , ecg indices %s "%(str(channel_name),self.ecg_indices))
            #self.ecg_epochs.plot_image(combine='mean')         
            self.all_ecg_scores[p] = self.ecg_scores
            self.all_ecg_indices[p] = self.ecg_indices      

        logger(inspect.stack())  
        print("ECG indices: ", self.all_ecg_indices)
        try:
            self.ecg_bads = []
            self.ecg_bad_ct = {}
            for index_all_channels in self.all_ecg_indices.values():
                for ica_index in index_all_channels:
                    self.ecg_bads.append(ica_index)
                    if ica_index in self.ecg_bad_ct.keys():
                        self.ecg_bad_ct[ica_index] = self.ecg_bad_ct[ica_index] +1
                    else:
                        self.ecg_bad_ct[ica_index] = 1
            self.ecg_bads = set(self.ecg_bads)
            self.ecg_bads = list(self.ecg_bads)
            print("ALL ECG bads ", self.ecg_bads)
            print("ECG BAD Counts: ", self.ecg_bad_ct)

            # Get two most frequently identified ICA indices to exclude
            max_key = max(self.ecg_bad_ct, key= self.ecg_bad_ct.get) if len( self.ecg_bad_ct) > 0 else None
            second_largest_key = None
            ica_counts = self.ecg_bad_ct.values()
            if len(self.ecg_bad_ct) > 1:
                second_largest = sorted(ica_counts)[max(-2,-len( self.ecg_bad_ct))]
            else:
                second_largest = None
            for count_num in ica_counts:
                for ica_index in self.ecg_bad_ct.keys():
                    if self.ecg_bad_ct[ica_index] == second_largest:
                        second_largest_key = ica_index
            print("Max ", max_key, " 2nd: ", second_largest_key)
            self.ica.excludes['ECG'] = [max_key,second_largest_key]

            # Exclude all identified components
            # self.ica.exclude = self.ica.exclude + self.ecg_indices        
            # Exclude the two most frequent ICA components, if there are enough channels identifying this component
            if second_largest_key != None and self.ecg_bad_ct[second_largest_key] >= min_channels_ica:
                self.ica.exclude = self.ica.exclude + [second_largest_key]
            if max_key!=None and self.ecg_bad_ct[max_key] >= min_channels_ica:
                self.ica.exclude = self.ica.exclude + [max_key]
        except:
            self.warn("Issue getting ECG Indices")  


        if self.ecg_bads:
            #self.ica.plot_sources(self.raw_cleaned, picks=self.ecg_indices)
            #self.ica.plot_properties(self.raw_cleaned, picks=self.ecg_indices, psd_args={'fmax': 35.})
            #self.ica.plot_sources(self.epochs_ransac_autoreject, picks=self.ecg_indices);
            #self.ica.plot_properties(self.epochs_ransac_autoreject, picks=self.ecg_indices, psd_args={'fmax': 35.});
            try:
                self.ecg_evoked = mne.preprocessing.create_ecg_epochs(self.raw_cleaned).average()
                #self.ecg_evoked.plot_image(combine='mean')
                #self.ecg_evoked.plot(gfp='only', spatial_colors=True, ylim=dict(eeg=[-12, 12]))
            except Exception as e:
                self.warn("Issue getting ECG Evoked " +str(e))  
            # try:
            #     mne.viz.plot_epochs(self.ecg_evoked)
            # except:
            #     self.warn("Cannot print ECG Evoked")  

            # self.ica.plot_scores(self.ecg_scores);

        try:
            # Plot ECG
            self.ecg_avg_df = self.ecg_evoked.to_data_frame()
            
            #self.ecg_avg_df['sum'] = self.ecg_avg_df.sum(axis=1) 
            #Rather, sum up excluding bad channels
            self.ecg_avg_df = sum_df(self.ecg_avg_df)

            # Plot ERFs
            # plt.plot(self.ecg_avg_df['time'].values, self.ecg_avg_df['sum'].values)
            # plt.legend(["ECG plot"], bbox_to_anchor=(0.75, 1.15), ncol=2)
            # plt.show()
        except Exception as e:
            self.warn("Cannot print ECG Epochs " +str(e))  


    def print_conditions(self):
        try:  
            print("CONDITIONS ARE: ", self.cond_A, self.cond_B)
            print("Conds we care about ", self.conds_we_care_about)
            print("Epochs post ransac autoreject", self.epochs_ransac_autoreject)
            print("Epochs generic ", self.evoked_generic)
        except:
            self.warn("No conditions defined")   
        

    def apply_ICA(self):
        '''Applies ICA fits'''

        #print("Plotting data before removal")
        #self.raw_cleaned.plot(start=START_TIME, duration=DURATION,n_channels=N_CHANS);
        if hasattr(self.ica, 'exclude'): # self.ica.exclude:
            logger(inspect.stack())  
            print("ICA components to remove: ", self.ica.exclude)
            self.ica.exclude = np.unique(self.ica.exclude).tolist()
            logger(inspect.stack())  
            self.ica.apply(self.raw_cleaned) # In-place replacement
            logger(inspect.stack())  
            self.ica.apply(self.epochs_ransac_autoreject) # In-place replacement

            # plot diagnostics
            if self.ica.exclude:
                pass
                #self.ica.plot_properties(self.raw_cleaned, picks=self.ica.exclude, psd_args={'fmax': 35.});
                #self.ica.plot_properties(self.epochs_ransac_autoreject, picks=self.ica.exclude, psd_args={'fmax': 35.});

                #self.cleaned_ica = self.raw_cleaned.interpolate_bads(reset_bads=True)
                #print("Plotting ICA after removal")      
                #self.ica.plot_sources(self.epochs_ransac_autoreject, show_scrollbars=False);        
                #self.ica.plot_sources(self.raw_cleaned, show_scrollbars=False);    
                #self.ica.plot_overlay(self.raw_cleaned, exclude=self.ica.exclude, picks='mag');    # ,show=False    

                #print("Plotting data after removal")
                #self.raw_cleaned.plot(start=START_TIME, duration=DURATION,n_channels=N_CHANS);
                
            else:
                print("Not excluding anything")

    def examine_ica(self):
        '''
        Run some plots to show what ICA has done
        '''
        try:
            self.describe()
            self.plot_ecg()
            self.plot_eog()

            print("Plotting before ICA")
            self.ica.plot_sources(self.raw_cleaned_copy, show_scrollbars=False);  
            self.raw_cleaned_copy.plot(start=120, duration=30,n_channels=20,remove_dc=False)
            plt.show()
            
            #self.ica.plot_sources(self.raw_cleaned, show_scrollbars=False); 
            self.ica.plot_sources(self.epochs_ransac_autoreject, show_scrollbars=False);             
            
            print("Plotting after ICA")
            self.epochs_ransac_autoreject.plot(start=120, duration=30,n_channels=20,remove_dc=False)
            #self.raw_cleaned.plot(start=120, duration=30,n_channels=20,remove_dc=False)
            plt.show()     
            
        except Exception as e:
            self.warn("ICA object missing " + str(e))
        

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
            if ch_name not in self.epochs_ransac_autoreject.info['ch_names']: #self.raw.info['ch_names']:
                self.missing_channels.append(ch_name)
            c+=1
        self.missing_channels = list(set(self.missing_channels))

    def remove_missing_channels(self):
        '''Remove channels from all participant-associated objects based on a ground-truth of which channels should be there'''

        self.identify_missing_channels()
        for missing_channel in self.missing_channels:
            status = missing_channel in self.epochs_ransac_autoreject.info['ch_names']
            if status:
                print("Deleting channel 2, ", missing_channel)
                for attrib in ['raw_cleaned','evoked_generic','evoked_all','epochs','epochs_ransac','epochs_ransac_autoreject']:
                    if hasattr(self,attrib):
                        obj = getattr(self, attrib)
                        obj.drop_channels([missing_channel]) 
            else:
                print(missing_channel + " was not found, not removing")
        self.signals_of_interest()

    def save_participant(self, cropped=True):
        '''Save participant object as pickle'''

        ## Determine name to save based on what was run
        string_save = " " + events_to_tag[1] + " "
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
        if cropped:
            string_save += " " +str(int(SEC_TO_MS*EPOCH_START_ANALYSIS)) + "to" +str(int(SEC_TO_MS*EPOCH_END_ANALYSIS))
        else:
            string_save += " " +str(int(SEC_TO_MS*EPOCH_START)) + "to" +str(int(SEC_TO_MS*EPOCH_END))


        if CONT_PREDICTOR_USE_DIFF_CUTOFF:
            string_save += " "+str(int(100/NUM_BINS_SURPRISE)) + "%"+ "vs"+str(int(100/NUM_BINS_SURPRISE)) + "%"
        else:
            string_save +=" " +str(int(100/NUM_BINS_SURPRISE)) + "%"+ "vs"+str(100-int(100/NUM_BINS_SURPRISE)) + "%"
        print("File to be saved as :", string_save)
        SAVE_EXPERIMENT_NAME = string_save

        wait_until_memory_free(required_memory = 3, max_wait_time_mins = 1.5) # This requires some memory (set conservatively here as I may run many threads, should be < 500MB)
        # Dump participant by participant
        pickle.dump( self, open( BASE_FOLDER+'Experiments Final\\' + str(self.p_id) +string_save+str(".pickle"), "wb" ) )
        
        
    def get_surprise_dir(self):
        '''
        Get the folder in which the predictors are stored
        '''
        if float(self.p_id) >= CUTOFF_FILE_NAMES:  # Same file for all participants after 9000.
            #self.surprise_file_dir = r"E:\BigData\\MEG\\MRES\\ME125_MMN_phase1_Yanan\\Common\\"
            self.surprise_file_dir = r"E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\HGF Settings Used\\"
        else:
            self.surprise_file_dir = self.BASE

    def get_surprise_hgf(self):
        '''
        Load the predictor of surprise for the HGF models
        '''
#         Loads saved expected surprise for each event
        
#         Vanilla    HGF     stored e.g. OLDER PTCPS E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\2689\\HGF346 PE2_new.txt
#         Vanilla    HGF     stored e.g. NEWER PTCPS E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Common\\HGF346 PE2_new.txt
        
#         Integrated HGF old stored e.g. OLDER PTCPS E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\2689\\HGF346 PE2_mod_baked_new.txt
#         Integrated HGF old stored e.g. NEWER PTCPS E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\HGF Settings Used\\HGF346 PE2_mod_baked_new.txt

#         Integrated HGF new stored e.g. OLDER PTCPS E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\2689\\HGF346 PE2_mod_baked v2_new.txt
#         Integrated HGF new stored e.g. NEWER PTCPS E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\HGF Settings Used\\HGF346 PE2_mod_baked v2_new.txt
        
        for time in HGF_TIMES:
            time = str(time)
            for hgf_type in ['PE2','PE3','PWPE2','PWPE3','PE2_mod','PE2_mod_baked','PE3_mod_baked','PWPE2_mod_baked','PWPE3_mod_baked']:
                
                if "mod_baked" in hgf_type:   # Integrated HGF new
                    directory = self.surprise_file_dir
                    full_path = directory+'HGF'+time + " " + hgf_type +" v2_new.txt"                    
                else:                   # Vanilla HGF
                    directory = self.surprise_file_dir
                    full_path = directory+'HGF'+time + " " + hgf_type +"_new.txt"

                #print("Hgf type %s , opening %s:"%(str(hgf_type),str(full_path)))
                try:               
                    with open(full_path) as f:  
                        if not hasattr(self,'surprise'):
                            self.surprise = {}
                        self.surprise['HGF_'+time+'_' + hgf_type] = f.readlines()
                except:
                    self.warn("Cannot find HGF "+str(time)+" " + hgf_type + " for: " + str(self.p_id) + " Looked at : " + str(full_path))   


    def get_surprise_BS(self):
        '''Loads saved expected Bayesian surprise for each event
        '''
        try:
            full_path = self.surprise_file_dir+ 'BS'+RELEVANT_LENGTH_PREDICTOR+'_new.txt'
            with open(full_path) as f:
                self.surprise['BS'] = f.readlines()
        except:
            self.warn("Cannot find Bayesian surprise file for: " + str(self.p_id) + " Looked at : " + str(full_path))

    def get_surprise_CS(self):
        '''Loads saved expected Confidence-Corrected Surprise for each event
        '''
        try:
            full_path = self.surprise_file_dir+ 'CS'+RELEVANT_LENGTH_PREDICTOR+'_new.txt'
            with open(full_path) as f:
                self.surprise['CS'] = f.readlines()
        except:
            self.warn("Cannot find Confidence-Corrected surprise file for: " + str(self.p_id) + " Looked at : " + str(full_path))     
    def get_surprise_PS(self):
        '''Loads saved expected Predictive surprise for each event
        '''
        try:
            full_path = self.surprise_file_dir+ 'PS'+RELEVANT_LENGTH_PREDICTOR+'_new.txt'
            with open(full_path) as f:
                self.surprise['PS'] = f.readlines()
        except:
            self.warn("Cannot find Predictive Surprise file for: " + str(self.p_id) + " Looked at : " + str(full_path))

    def get_surprise(self):
        '''
        Loads the expected surprise using various models for the participant
        '''
        self.surprise = {}
        self.get_surprise_dir() # Run before to know the directory to look in for the participant
        self.get_surprise_hgf()
        self.get_surprise_BS()
        self.get_surprise_CS()      
        self.get_surprise_PS()
        
        for key in self.surprise.keys():
            self.surprise[key] = [float(item.replace("\n","")) for item in self.surprise[key]]
        self.get_surprise_distributions()
        
    def get_surprise_distributions(self,equal_width=False, NUM_BINS=NUM_BINS_SURPRISE):
        '''Calculates a distribution of the values a predictor takes
        Then spits out NUM_BINS different predictor values, with each bin containing the same number of samples

        NOTE the effect of NUM_BINS
        high surprise cutoff is n-2th item out of n bins
        low surprise cutoff is 2nd item out of n bins
        e.g.
        NUM_BINS = 5 ->  high surprise is 4th out of 5, i.e top 20%,    low surprise is 2nd out of 5, i.e bottom 20%
        NUM_BINS = 10 -> high surprise is  9 out of 10, i.e. top 10%),  low surprise is 1st out of 10, i.e. bottom 10%        
        '''
        
        self.surprise_cutoffs = {}
        self.surprise_cutoffs_low = {}
        for surprise_key in self.surprise.keys():
            if equal_width:
                min_ = float(min(self.surprise[surprise_key]))
                max_ = float(max(self.surprise[surprise_key]))
                print("Surprise min and max: "+str('{0:.2f}'.format(min_),'{0:.2f}'.format(max_)))
                intervals = np.arange(min_,max_,(max_-min_)*0.1)
                j = 0
                cts = {}
                for val in self.surprise[surprise_key]:
                    for interval_start in intervals:
                        interval_end = intervals[min(j+1,len(intervals)-1)]
                        val = float(val)
                        if val > interval_start and val <= interval_end:
                            if '{0:.2f}'.format(interval_start) in cts:
                                cts['{0:.2f}'.format(interval_start)] = cts['{0:.2f}'.format(interval_start)]+1
                            else:
                                cts['{0:.2f}'.format(interval_start)] = 1 
                    j+=1
                self.surprise_cutoffs[surprise_key] = find_cutoff(cts)
                print("Finding cutoff...")
            else: # Then use equal frequencies
                data = self.surprise[surprise_key]

                #create histogram with equal-frequency bins 
                n, bins, patches = plt.hist(data, equalObs(data, NUM_BINS), edgecolor='black');
                #plt.show()

                #display bin boundaries and frequency per bin 
                #print(bins, n)
                cutoff = bins[len(bins)-2]  # high surprise cutoff (n-2 item out of n bins e.g. 4th out of 5, i.e top 20%,    or  9 out of 10, i.e. top 10%)
                cutoff_low = bins[1]        # low surprise cutoff (2nd item out of n bins e.g. 2nd out of 5, i.e bottom 20%,    or  1 out of 10, i.e. bottom 10%)
                
                #print("Surprise cutoffs, ", surprise_key, cutoff, cutoff_low)
                self.surprise_cutoffs[surprise_key] = cutoff
                self.surprise_cutoffs_low[surprise_key] = cutoff_low

            
        #print("Surp cutoffs: ", self.surprise_cutoffs)    
        #print("Surp cutoffs low: ", self.surprise_cutoffs_low) 
    

    def apply_z_score(self,STANDARDISE=True,channel_wise=False):
        if STANDARDISE:
            for attrib in ['raw_cleaned','evoked_generic','evoked_all','epochs','epochs_ransac','epochs_ransac_autoreject']:
                if hasattr(self,attrib):
                    obj = getattr(self, attrib)
                    if channel_wise == True:
                        if type(obj) == type({}):
                            keys = list(obj.keys())
                            for key in keys:
                                print("Applying, ", key)
                                obj[key].apply_function(fun=calc_z_score, n_jobs=4)
                        else:
                            obj.apply_function(fun=calc_z_score, n_jobs=4)
                    else:
                        try:
                            obj.apply_function(fun=calc_z_score, n_jobs=4, channel_wise=False)  
                        except Exception as e:
                            keys = list(obj.keys())
                            for key in keys:
                                print("Applying, ", key)
            self.signals_of_interest()
            self.applied_z_score = True
        else:
            self.applied_z_score = False
    



    def create_evoked_dfs(self):

        '''Create evoked DFs'''
        # For all
        self.evoked_all = self.epochs_ransac_autoreject.average()
        self.evoked_all_df = self.evoked_all.to_data_frame()
        # For nth rep
        self.evoked_generic = {}
        self.evoked_generic_dfs = {}

        cond_key = self.conds_we_care_about.keys()
        for j in cond_key:       
            self.evoked_generic[str(j)] = self.epochs_ransac_autoreject[str(j)].average() 
            # Create df of evoked responses for the epoch
            self.evoked_generic_dfs[str(j)] = self.evoked_generic[str(j)].to_data_frame()

    def plot_evoked_dfs(self):
        '''
        Plot all conditions' ERF
        '''

        #self.evoked_all_df['sum'] = self.evoked_all_df.sum(axis=1)
        #Rather, sum up excluding bad channels
        self.evoked_all_df = sum_df(self.evoked_all_df)#,bads=self.bads)

        # plt.plot(self.evoked_all_df['time'].values, self.evoked_all_df['sum'].values)
        # plt.legend(["Average ERF for ALL conditions"], bbox_to_anchor=(0.75, 1.15), ncol=2)
        # plt.show()
        str_keys = self.conds_we_care_about.keys()
        for j in str_keys:       

            #self.evoked_generic_dfs[j]['sum'] = self.evoked_generic_dfs[j].sum(axis=1)
            self.evoked_generic_dfs[j] = sum_df(self.evoked_generic_dfs[j])#, bads=self.bads)

            # plt.plot(self.evoked_generic_dfs[j]['time'].values, self.evoked_generic_dfs[j]['sum'].values)
            # plt.legend(["Average ERF for condition "+j], bbox_to_anchor=(0.75, 1.15), ncol=2)
            # plt.show()

    def create_mmr_dfs(self):
        '''Create MMR DFs of an epoch for of a specific condition vs *ALL*
        '''

        plot_legend = []
        # For nth rep        
        self.mmr_dfs = {}
        str_keys = self.conds_we_care_about.keys()
        for j in str_keys:
            #print("CONDITION IS: ",str(j))
            #print("Event tag is :", self.conds_we_care_about[j])
            self.mmr_dfs[j] = copy.deepcopy(self.evoked_generic_dfs[j])
            # Calculates the MMR
            for col in self.mmr_dfs[j].columns.values:
                if col in self.evoked_generic_dfs[j].columns.values and col in self.evoked_all_df.columns.values:
                    self.mmr_dfs[j][col] = self.evoked_generic_dfs[j][col] - self.evoked_all_df[col] # When there are few conditions which appear equally often, these lines will be negatively correlated across various js
                else:
                    print("Column not in one of the dfs 1...")

            # Set time to the correct value manually  since we dont take a mismatch of it!
            self.mmr_dfs[j]['time'] = self.evoked_generic_dfs[j]['time']

            # Plot MMR over all channels
            #self.mmr_dfs[j]["sum"] = self.mmr_dfs[j].sum(axis=1)
            #Rather, sum up excluding bad channels
            self.mmr_dfs[j]["sum"] = sum_df(self.mmr_dfs[j])#, bads=self.bads)


    def create_mmr_difference_dfs(self,n_1,n_2):
        '''Create and plot ERF MMR'''

        if str(n_1) in self.conds_we_care_about.keys() and str(n_2) in self.conds_we_care_about.keys():
            # n1_th and n2_th deviant, plot both
            plot_legend = []

            # plt.plot(self.mmr_dfs[str(n_1)]['time'].values, self.mmr_dfs[str(n_1)]['sum'].values)
            # plot_legend.append("MMR after "+str(n_1)+ " condition")
            # plt.plot(self.mmr_dfs[str(n_2)]['time'].values, self.mmr_dfs[str(n_2)]['sum'].values)
            # plot_legend.append("MMR after "+str(n_2)+ " condition")      
            # plt.legend(plot_legend, bbox_to_anchor=(0.75, 1.15), ncol=2)
            # plt.show()        

            self.difference_df = self.evoked_generic_dfs[str(n_1)].copy()
            for col in self.mmr_dfs[str(n_1)].columns.values:
                if col not in ['time',"sum"]:
                    if col in self.evoked_generic_dfs[str(n_1)].columns.values and col in self.evoked_generic_dfs[str(n_2)].columns.values:
                        self.difference_df[col] = self.evoked_generic_dfs[str(n_1)][col] - self.evoked_generic_dfs[str(n_2)][col]
                    else:
                        print("Column not in one of the dfs 2...")

            
            ##Plot MMR per channel
            # for col in self.difference_df.columns.values[1:10]:
            #     print(col)
            #     plt.plot(self.difference_df[col]['time'].values, self.difference_df[col].values) 


            # MMR over all channels
            # self.difference_df["sum"] = self.difference_df.sum(axis=1) # Includes bad channels and irrelevant columns (e.g. non MEG channels, the 'time' column etc)
            # Rather, sum up excluding bad channels
            self.difference_df = sum_df(self.difference_df)#, bads=self.bads)
            # plt.plot(self.difference_df['time'].values, self.difference_df['sum'].values)    
            # plt.legend([str(n_1) +" MMR minus " + str(n_2)+" MMR"], bbox_to_anchor=(0.75, 1.15), ncol=2)
            # plt.show()
            
        else:
            print("One of the two passed arguments, ", n_1, " or ", n_2, " are not in the conditions dictionary:" )
            print("Conditions dictionary :", self.conds_we_care_about.keys())


    def epochs_trigger_leakage(self):

        '''Epoching for examining trigger leak'''
        self.find_events()
        
        self.events_processing_define_events()

        self.epochs_trigger_leakage = mne.Epochs(
                                                self.raw_cleaned,
                                                self.events,
                                                event_id=self.conds_we_care_about,
                                                tmin=-0.10, # Peri-stimulus period to check for trigger leakage
                                                tmax= 0.10,
                                                baseline=(None,None),
                                                #picks=None,
                                                #preload=True
                                                )

        ## DON'T Resample when examining trigger leakage
        # self.epochs_trigger_leakage = self.epochs.resample(DOWNSAMPLE_FREQ)

        ## Choose everything - unecessary here but can be used to filter classes
        self.picks = mne.pick_types(self.epochs.info, meg="mag", eeg=False, stim=False, eog=False, include=[], exclude=[])


    def evoked_trigger_leakage(self):
        '''Average ERFs around the trigger time
        '''
        self.evoked_trigger_leak = self.epochs_trigger_leakage.average() 

    def process_trigger_leakage(self):
        '''Produce an Evoked that is the average ERF around the time of trigger'''
        self.find_events() # We do not adjust the event times for the sound, as we are just trying to see if the MEG trigger casues leakage into the MEG channels
        self.epochs_trigger_leakage()
        self.evoked_trigger_leakage()
        
    def match_channels(self):
        '''Matches the adult sensors to appropriate child sensors 
        So we leave the child channels as 1,2,3,...125 and then we pick the adult channels that best match these
        This is necessary if e.g. plotting the topography of the brain signal as MEG001 in the adult system does not correspond to MEG001 in the child system (for example) in terms of where it is located in the brain

        It's not enough to just drop the non-comparable channels; we also have to fool MNE into thinking that the channels are comparable when e.g. doing a grand average as otherwise it'll drop the channels
        To do this we need to give the corresponding adult channels the same names as the children's (i.e. drop 126 to 160)
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
      
    def remove_surplus_channels(self):
        '''Search and delete channel from any relevant object'''
        required_channels = []
        for ch in range(1,126): # MEG 001 to MEG 125 (MEG 071 will be missing in both child and adult systems, so 124 channels will be left over after this)
            ch_name = "MEG "+str(ch).zfill(3)
            required_channels.append(ch_name)
        self.surplus_channels = []
        ptcp_chs = self.epochs_ransac_autoreject.info['ch_names']
        for ptcp_ch in ptcp_chs:
            if ptcp_ch not in required_channels and "MEG" in ptcp_ch:
                print("Ptcp ID %s chan %s is surplus, dropping"%(str(self.p_id), str(ptcp_ch)))
                self.surplus_channels.append(ptcp_ch)

        for attrib in ['raw_cleaned','evoked_generic','evoked_all','epochs','epochs_ransac','epochs_ransac_autoreject']:
            if hasattr(self,attrib):
                obj = getattr(self, attrib)
                if type(obj) == type({}):
                    keys = list(obj.keys())
                    for key in keys:
                        for chan in self.surplus_channels:
                            status = chan in obj[key].info['ch_names']
                            if status and "MEG" in chan:
                                print("Deleting channel 1 obj key, ", chan)
                                obj[key].drop_channels([chan]) 
                            else:
                                print(chan + " was not found, not removing")
                else:
                    for chan in self.surplus_channels:
                        status = chan in obj.info['ch_names']
                        if status  and "MEG" in chan:
                            print("Deleting channel 1 obj, ", chan)
                            obj.drop_channels([chan]) 
                        else:
                            print(chan + " was not found, not removing")
        self.signals_of_interest()

    def pick_and_rename(self,obj,chans_to_keep,mapping,rename=True):
        obj.pick_channels(chans_to_keep)
        if rename:
            obj.rename_channels(mapping)

    def align_channels_paul_individual(self):
        '''Perform a remapping of the adult channel names to child channel names to allow comparison across different machines
        This function doesn't do the deletion, this happens elsewhere
        '''
        

        logger(inspect.stack()) 
        if self.is_adult_system:
            print("H1", len(self.epochs_ransac_autoreject.info['ch_names']))
            print("Scanned on adult system, so remapping channels for ptcp: ", self.p_id)
            data_dir = r"E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Lance\MEG\Python MEG\Change_channel_mapping\\"
            MEG125_raw = mne.io.read_raw_kit(data_dir + "MEG125_sample.con",
                                                allow_unknown_format=False,
                                                verbose=True)
            #print(MEG125_raw.info)
            
            data_dir = r"E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Lance\MEG\Python MEG\Change_channel_mapping\\"
            # load up the mapping lists - stored as .mat files
            MEG125 = scipy.io.loadmat(data_dir + "labels_systemmatch_child_sensors.mat", simplify_cells=True)
            MEG160 = scipy.io.loadmat(data_dir + "labels_systemmatch_adult_sensors.mat", simplify_cells=True)    
            
            
            MEG125_picks = MEG125["lab125"]
            MEG160_picks = MEG160["lab160"]

            for idx, ch in enumerate(MEG160_picks):
                MEG160_picks[idx] = ch.replace("AG", "MEG ")
            for idx, ch in enumerate(MEG125_picks):
                MEG125_picks[idx] = ch.replace("AG", "MEG ") 

            MEG125_picks = list(set(MEG125_picks)-set(['MEG 072']))
            #MEG160_picks = list(set(MEG160_picks)-set(['MEG 072']))
            print("H2", len(MEG125_picks))
            print("H3", len(MEG160_picks))    
            # for m in range(126,161):
            #     MEG160_picks = list(set(MEG160_picks)-set(["MEG "+str(m).zfill(3)]))
            print("H4", len(MEG160_picks))
            
            MEG125_raw.pick_channels(MEG125_picks)
            print("H5a", len(MEG125_raw.info['ch_names']))


            self.get_raw(trigger_channels_known=False)
            ## Finding events
            events = mne.find_events(
                self.raw,
                output="onset",
                consecutive=False,
                min_duration=0,
                shortest_event=1,  # 5 for adults
                mask=None,
                uint_cast=False,
                mask_type="and",
                initial_event=False,
                verbose=None,
            )


            #MEG160_epochs = mne.Epochs(MEG160_raw, events, tmin=-0.1, tmax=0.4, preload=True)
            #MEG160_epochs.pick_types("mag")

            print("H5b", len(sorted(self.epochs_ransac_autoreject.info['ch_names'])))

            #info = mne.create_info(ch_names, self.epochs_ransac_autoreject.info['sfreq'], ch_types='mag', verbose=True)
            for attrib in ['raw_cleaned','evoked_generic','evoked_all','epochs','epochs_ransac','epochs_ransac_autoreject']:
                if hasattr(self, attrib):
                    obj = getattr(self, attrib)
                    if type(obj) == type({}):
                        keys = list(obj.keys())
                        for key in keys:

                            try:
                                print("H6", len(obj[key].info['ch_names']) )
                                obj[key].pick_types("mag")
                                print("H7", len(obj[key].info['ch_names']) )
                            except:
                                print("PROBLEM PICKING MAG A", attrib, key)
                            try:
                                obj[key].pick_channels(MEG160_picks)
                                print(self.p_id, self.is_adult_folder, self.is_adult, self.epochs_ransac_autoreject)
                                print("H8", len(obj[key].info['ch_names']) )
                            except:
                                print("PROBLEM PICKING CHANNELS A", attrib, key)
                                print("H9", "MEG 160 picks ", len(MEG160_picks), sorted(MEG160_picks))
                                for ch in MEG160_picks:
                                    if ch not in obj[key].info['ch_names']:
                                        print("Mismatch 1", ch)
                                for ch in obj[key].info['ch_names']:                                    
                                    if ch not in MEG160_picks:
                                        print("Mismatch 2", ch)
                            try:         
                                obj[key].reorder_channels(MEG160_picks)
                            except:
                                print("PROBLEM REORDERING CHANNELS A", attrib, key)
                                print("MEG 160 picks ", len(sorted(MEG160_picks)))
                            try:
                                obj[key].info["chs"] = MEG125_raw.info["chs"]
                            except:
                                print("PROBLEM SETTING CHS A", attrib, key)
                                print(sorted(obj[key].info["chs"]))
                                print(sorted(MEG125_raw.info["chs"]))
                            try:
                                obj[key].info["ch_names"] = MEG125_raw.info["ch_names"]
                            except:
                                print("PROBLEM SETTING CH NAMES A", attrib, key)
                                print(sorted(obj[key].info["chs"]))
                                print(sorted(MEG125_raw.info["chs"]))

                    else:
                        # obj.pick_types("mag")
                        # obj.pick_channels(MEG160_picks)
                        # obj.reorder_channels(MEG160_picks)
                        # obj.info["chs"] = MEG125_raw.info["chs"]
                        # obj.info["ch_names"] = MEG125_raw.info["ch_names"]
                        try:
                            print("H6b", len(obj.info['ch_names']), sorted(obj.info['ch_names']))
                            obj.pick_types("mag")
                            print("H7b", len(obj.info['ch_names']), sorted(obj.info['ch_names']))
                        except:
                            print("PROBLEM PICKING MAG B", attrib)
                        try:
                            obj.pick_channels(MEG160_picks)
                            print("H8b", len(obj.info['ch_names']), sorted(obj.info['ch_names']))
                        except:
                            print("PROBLEM PICKING CHANNELS B", attrib)
                            print("H9b", "MEG 160 picks ", len(MEG160_picks), sorted(MEG160_picks))
                            for ch in MEG160_picks:
                                if ch not in obj.info['ch_names']:
                                    print("Mismatch 1", ch)
                            for ch in obj.info['ch_names']:                                    
                                if ch not in MEG160_picks:
                                    print("Mismatch 2", ch)
                        try:           
                            print("H10b ", sorted(MEG160_picks))
                            obj.reorder_channels(MEG160_picks)
                            print("H11b ", sorted(MEG160_picks))
                        except:
                            print("PROBLEM REORDERING CHANNELS B", attrib)
                            print("MEG160 picks ", sorted(MEG160_picks))
                            print(sorted(self.epochs_ransac_autoreject.info['ch_names']))
                        try:
                            obj.info["chs"] = MEG125_raw.info["chs"]
                        except:
                            print("PROBLEM SETTING CHS B", attrib)
                            print(sorted(obj.info["chs"]))
                            print(sorted(MEG125_raw.info["chs"]))
                        try:
                            obj.info["ch_names"] = MEG125_raw.info["ch_names"]
                        except:
                            print("PROBLEM SETTING CH NAMES B", attrib)
                            print(sorted(obj.info["chs"]))
                            print(sorted(MEG125_raw.info["chs"]))

            try_del_attr(self, 'raw')

    def align_channels(self,mapping=None):
        '''Perform a remapping of the adult channel names to child channel names to allow comparison across different machines
        This function doesn't do the deletion, this happens elsewhere
        '''

        if mapping==None:
            self.match_channels()
            mapping = self.sensor_mapping_dct_simple

        adult_surviving_chans = list(mapping.keys())
        child_surviving_chans = list(mapping.values())
        
        # logger(inspect.stack())
        # print(hasattr(self,'epochs'))
        # if hasattr(self,'epochs'):
        #     print(dir(self.epochs))
        #ch_names = [x for x in self.raw_cleaned.info.ch_names if 'MEG' in x]
        ch_names = [x for x in self.epochs_ransac_autoreject.info['ch_names'] if 'MEG' in x]
        #ch_names = [x for x in self.epochs.info['ch_names'] if 'MEG' in x] 
        #logger(inspect.stack())

        #info = mne.create_info(ch_names, self.epochs_ransac_autoreject.info['sfreq'], ch_types='mag', verbose=True)
        for attrib in ['raw_cleaned','evoked_generic','evoked_all','epochs','epochs_ransac','epochs_ransac_autoreject']:
            if hasattr(self, attrib):
                obj = getattr(self, attrib)
                if type(obj) == type({}):
                    keys = list(obj.keys())
                    for key in keys:
                        if len(ch_names) == num_meg_chs_adult:
                            self.pick_and_rename(obj=obj[key],chans_to_keep=adult_surviving_chans,mapping=mapping,rename=True)
                        elif len(ch_names) == num_meg_chs_child:
                            self.pick_and_rename(obj=obj[key],chans_to_keep=child_surviving_chans,mapping=mapping,rename=False)
                else:
                    if len(ch_names) == num_meg_chs_adult:
                        self.pick_and_rename(obj=obj,chans_to_keep=adult_surviving_chans,mapping=mapping,rename=True)
                    elif len(ch_names) == num_meg_chs_child:
                        self.pick_and_rename(obj=obj,chans_to_keep=child_surviving_chans,mapping=mapping,rename=False)

        self.remove_surplus_channels()
        # These should all agree (make it an assert?)
    #     print(len(self.evoked_generic[condition_A].info['ch_names']))
    #     print(len(self.evoked_all.info['ch_names']))    
    #     print(len(self.epochs_ransac_autoreject[condition_A].info['ch_names']))        
