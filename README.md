To run:  

- i) Install packages as per requirements.txt  
- ii) Set **config.py**  
- iii) Update bad participants and channels per participant ID **(participant_data.py**)  
- iii) Update directories as needed in main class objects (**experiment.py, participant.py**) and in **standardise_all_3.py**    
- <Note: I'm going to replace these steps with a **do_everything.py ** >  
-- iv) Run **process_all_participants.py**  
--- Note that this can be sped up by using multiple copies (on different cores) and changing file_num_start and file_num_end to be mutually exclusive between the copies  -- v) Create a standardised copy of their data by running **standardise_all_3.py**  
-- vi) Run **save_output_graphs.py**   


Experiment code to run in the lab is in the **RunExperiment** folder  
HGF output predictions for my participants with IDs >= 9000 in the **Predictions** folder  
Other helpful files in the **Helpers** folder  
  
It is recommended to run on a machine with >=4GB of memory.  
24GB of memory was required when loading, analysing and saving the output for all 72 participants at once (1000Hz files).  
If this is not available, you can run steps (i-iv), and then (v, vi) separately once (i-iv) is concluded  

