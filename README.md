To run (Note: I'm going to replace steps 5-7 with a **do_everything.py **):    

1) Install packages as per requirements.txt  
2) Set **config.py**  
3) Update bad participants and channels per participant ID **(participant_data.py**)  
4) Update directories as needed in main class objects (**experiment.py, participant.py**) and in **standardise_all_3.py**    
5) Run **process_all_participants.py**, which saves each participant's processed data to disk  
- Note that this can be sped up by running multiple copies of process_all_participants (on different cores)  
- To do this, change file_num_start and file_num_end to be mutually exclusive between the copies  (e.g. (0, 12) on one, (12,24) on a 2nd core and (24,36) on a third  
6) Create a standardised copy of the saved participant data by running **standardise_all_3.py**  
7) Set parameters in and then run **save_all_graphs.py**   


  
  
  
  
Experiment code to run in the lab is in the **RunExperiment** folder (requires Cogent2000)
HGF output predictions for my participants with IDs >= 9000 in the **Predictions** folder  
Other helpful files are found in the **Helpers** folder  
  
It is recommended to run on a machine with >=4GB of memory.  
24GB of memory was required when loading, analysing and saving the output for all 72 participants at once (1000Hz files).  
If this is not available, you can run steps (1-5), and then (6-7) separately once (1-5) is concluded  

