
# Sets up participants

# Channels noted by the experimenter to be bad (see automated readout of *NOTES*.txt file) and altered upon inspection
# LANCE: in python had to reference 1 different to the notes file

#import helpers
from helpers import *
from config import *




CHILD_BASE          = BASE_FOLDER+'\Child_MEG\\'
ADULT_BASE          = BASE_FOLDER+'\Adult_MEG\\'


MANUAL_BAD = {   # HFN                                                  # Digital Noise
                '2683':        [],
                '2689':        list(range(80,90))                                   + [13,15,         31] + [42], 
                '2717':        list(range(82,86))                      + [51,53,54,55,56,57,159],
                '2718':        list(range(82,86))                      + [51,53,54,55,56,57,159],
                '2719': [24,            80,84,               88]                    + [119],
                '2722':                                             [94], # Trapped flux in many channels
                '2723':           [82] + list(range(84,90)),
                '2729':           [82] + list(range(84,90))                                          +[31],
                '2730':           [82] + list(range(84,90))                       + [119],
                '2733': [80,82] + list(range(84,90))         + [108],
                '2737': [24]           + list(range(81,90))            + [108,122,124],
                '2744': [24]           + list(range(81,90))            + [108],
                '2745': [24]           + list(range(81,90))            + [108]     + [13,14,15,16]   +[31,32],
                '2748': [24]           + list(range(81,90))  + [101],
                '2750': [24]           + list(range(81,90))  + [101],
                '2752': [24]           + list(range(81,90)),  
                '2760':                  list(range(80,90)) ,
              

                # HFN                          + Flatlined
                '3370': [45]                   +         [58],             
                '3374': [13,42,90]             + [20, 48, 58] +   [4,5,6,8,15,16.20,27,29,45,58,113,            122,123], # Last: Lance
                '3380':                                  [58],
                '3402': [20, 27, 41, 73, 87]   + [17, 20, 24, 25, 26, 58, 85, 86, 100] + list(range(88,93)),
                '3416': [24]                   + [20,     58],
                '3418': [25, 29]               +          [58,      62],
                '3420':      [44, 115]         + [20,     58],
                '3422':                          [20,     58],
                '3423':                                  [58],
                '3426': [88]                   + [20,     58]      + [13, 17, 20, 24, 34, 39, 42, 43] + [74,86], #  Last two: Saturations, Lance 
                '3429': [40, 118]              + [20,     58]      + [15, 16, 45],
                '3433':                                   [58],
                '3438':                                   [58],
                '3443':                          [20,      58,      62],
                '3448': [15, 45],
                '3508':     [45],
                '3612':     [45],
             }

MANUAL_BAD_LANCE = {'2552':                                                 list(range(81,90)) + [85],
                    '2589':  [29, 41] + [121, 122, 124, 129, 130, 152]      + list(range(81,90)),  
                    '2629':  [7, 11, 12, 40,    41,     49, 50 , 54, 55, 59, 60, 61, 63, 64, 66, 86, 100, 106,      109, 114,115,           122], #122 37Hz noise]]
                    '2632':  [3,7,11,12,13,16,21,25,26,32,34,35,40,41, 46,47,49,50,51, 53,54,55, 62, 64, 77, 100, 106, 109, 110, 115, 120,  122, 123],
                    '2678':  list(range(17,21)) + [24, 79]                  + list(range(81,90)),
                    '2683':  [109],
                    '2687':  [21, 41, 46, 49,     50, 102, 109, 115, 122, 123]  + list(range(81,90)),
                    '2689':  [17, 20, 24] + [40, 42]  # (Pulsey)
                                                                            + list(range(81,90)),
                    '2695':  [21,                   46, 49]                 + list(range(81,90)) +[102, 103,        109, 115]             +[122, 123],
                    '2696':  [2,                    46, 49,                             77, 78,                105, 109],
                    '2697':  [12, 16,               46, 49,     50,                             86, 89, 91, 115],
                    '2699':  [86, 109],
                    '2702':  [          33,     41, 46,                   52,                                       109],
                    '2703':  [2,        33,     41, 46, 49,       51, 52, 54, 55,     76,       88, 94, 105,        109],
                    '2712':  [21,       33, 39, 41, 46,           51 ,52,         75, 76,                           109],
                    '2713':  [      30, 33, 37, 41, 46, 47, 48,   51, 52],
                    '2716':  [21, 30, 37,                         51, 52,                                           109],     
                    '2717':  [7,15, 24]                                     + list(range(81,90)),
                    '2718':  [51,56,57] # Blocky,
                            + [152,153,156]                                 + list(range(81,90)),
                    '2719':  [24]                                           + list(range(81,90)) + list(range(97,104)) + [120,140],
                    '2723':  [20,       24,                     50,   52]   + list(range(81,90)) + list(range(97,103)) + [108,120,148],
                    '2729':  [          24]                                 + list(range(81,90)) + list(range(97,103)),
                    '2730':  [          24,                             79,                 81, 85, 89],


                    '2733':  [          24,                           52],
                    '2737':  [12,       24, 27, 29,                                         81, 85, 89, 90,                                 122,124],
                    '2738':   [     21,              46,                            75,                     91,      109],
                    '2744':  [20,       24, 27, 29,                                     79, 81, 85, 89, 90],
                    '2745':  [20,       24,                                                 81, 85, 89,                                     122,124],
                    '2748':  [20,       24,                                             79, 81, 85, 89,         95],
                    '2750':  [          24,                                                 81, 85, 89],
                    '2752':  [          24,                                                 81, 85, 89],   
                    '2760':  [20,       24,                                             79, 81, 85, 89],      
                    '2766':  [  21,                                                                                  109,                   122, 123],
                    '2785':  [                       46, 49,                55,                                      109, 110],
                    '2786':  [  21,                  46,                                                90,          109,                   122, 123],
                    '2787':  [16, 20, 22, 27,   41, 46, 49,     50,         55, 61, 63,     82, 86, 89, 91, 102, 105,109,   115, 120,       122, 123],
                    '2793':  [2,                41, 46, 49, 58, 77, 78, 89, 102, 109],
                    '2854':  [15,                   46,                                                              109],
                    '2858':  [21,                   46, 50, 109, 122, 123],  
                    '2866':  [                          49, 55,                                 87,                  109,                   122, 123],  
                    '2872':  [16, 17,               46, 49, 50, 55, 62, 63, 77, 78,                 89,    102, 105, 109,                   122, 123],
                    '2875':  [                      46,                                                              109,                   122, 123],  
                    '2888':  [                      46, 49, 54, 55, 62, 64, 67, 77, 78, 102, 105, 109, 123],
                    '2897':  [                      46,         55,                                              91, 109,                   122, 123],
                    '2908':  [                      46,         55,                                                  109,                   122, 123],
                    '2913':  [21, 22,               46, 109, 122, 123],  
                    '3380':  [15, 21, 33, 58 , 76],
                    '3394':  [58],
                    '3402':  [7, 33, 50, 52, 58, 63, 66, 74, 86, 89, 91, 105, 122, 123],
                    '3416':  [6, 58],
                    '3418':  [15, 20, 58],
                    '3419':  [20, 58, 116],
                    '3421':  [58, 74, 116, 119],
                    '3422':  [13, 14, 20, 58, 116],
                    '3423':  [15,58],
                    '3426':  [20, 58, 63, 74, 86, 88],
                    '3429':  [15, 45, 58, 71, 72, 86],
                    '3433':  [20, 22, 58],
                    '3434':  [22, 26, 66, 67, 79, 85],
                    '3438':  [58, 63, 66, 67],
                    '3439':  [15, 20, 58],
                    '3441':  [58],
                    '3443':  [33,46,50],
                    '3448':  [33, 50, 122, 123],
                    '3612':  [81, 86, 87, 92, 96],
                    '9000':  [12, 21, 22, 26, 27, 29, 50, 52, 58, 60,  67, 79, 95, 117, 139, 142, 144, 149, 154],
                    '9001':  [45, 46, 47, 48, 61, 62],
                  # Lance HFN
                    '9004': [7,10],
                    '9005': list(range(22,29)) + [94,117],
                    '9006': [11,                100],
                    '9008': [11, 40, 49,        100],
                    '9009': [11, 49, 89, 94,    100, 123, 148, 151, 153]     +  [21,23,    24,25,26,27,28], # Digital noise
                    '9010': [11, 21, 49, 89, 94, 100, 109, 123, 147, 151]           + [22, 24,  26,    28],
                    '9011': [11,     49, 89, 94, 100 ]                                      + [33],
                    '9012': [11, 21, 22, 26, 28, 49, 65, 66, 78, 89]                            + [71,77],  # Digital noise
                    '9013': [11, 20, 21, 44, 49] +list(range(80,95))  + [100, 117, 153],
                    '9017': [11, 20, 21,     49] + [64,65,66] + [78, 94, 100, 123, 153],
                    '9018': [11, 20, 21,     49] + [64,65,66] + [78, 94, 100, 123, 153],
                    '9019': [11,     21, 28, 49,                 89, 94,      123, 151],
                    '9020': [11,             49,                 89, 94, 100, 109, 141, 151, 153],
                    '9021': [6, 11,      40, 49,                         100,                     ]
                    }


# Pad out bad channels: replace with uniform list of bad channels, assuming these were not noted carefully for each participant
for key in range(2689,2760):
    if str(key) in MANUAL_BAD.keys():
        MANUAL_BAD[str(key)] = list(set(MANUAL_BAD[str(key)] + [24] +list(range(80,90)) + [108]))
for key in range(3370,3616):
    if str(key) in MANUAL_BAD.keys():
        MANUAL_BAD[str(key)] = list(set(MANUAL_BAD[str(key)] + [45] + [20, 58]))
                  

# Data manipulation
def try_add(dct,key,value):
    
    if key in dct.keys():
        #print(key,dct[key])
        dct[key] = dct[key] + [value]
    else:
        dct[key] = [value]   
for p_id in ['9000'] + list(range(9002,9010)): # Lance: All so far that I have done have had the same bad channels
    for chan in [11,21,49,89,100,109,123]:
        try_add(MANUAL_BAD, str(p_id), chan)
                        
# Merge dictionaries
for key in MANUAL_BAD_LANCE.keys():
    if key in MANUAL_BAD:
        MANUAL_BAD[key] = MANUAL_BAD[key] + MANUAL_BAD_LANCE[key]
    else:
        MANUAL_BAD[key] = MANUAL_BAD_LANCE[key]

# Remove duplicates   
for key in MANUAL_BAD.keys():
    MANUAL_BAD[key] = sorted(list(set(MANUAL_BAD[key])))       

# Give the exclusion channels names                      
for key in MANUAL_BAD:
    MANUAL_BAD[key] = ["MEG " + str(x).zfill(3) for x in  MANUAL_BAD[key]]








# Participants that were excluded from Hannah's experiment (ID < 9000) and mine (>9000)
# These have already been moved out of the folder, but explicitly exclude them in case of accidental re-inclusion in the folder
PROBLEM_PARTICIPANTS = ['2480',                     # .CON FILE ISN'T ~15 MINS (suspicious of what was shown, can check to add these back)  
                        '2607',                     # Error loading RAW, also there are up to 11 repetitions here, not 7
                        
                        '2640','2656','2701',       # Bad coreg 
                        '2751','3370','3412',       # Bad coreg 
                        '3440','3445','3505',       # Bad coreg 

                        '2681',                     # Huge noise on many channels e.g. on list(range(45,50)) + [86,93,101],
                        '2722',                     # Trapped flux, a lot of HF and LF noise
                        '2724','2810',              # No events on any sound channel that I can see.

                        #'2552','2689',             # No sound detected on Matlab pipeline -> found on other channels.
                        #'2733','2748' '2719',      # Lag between event and sound is extremely variable (especially at beginning) and had to be replaced all by medians
                        '3374',                     # Weird U-shaped bands found frequently on many channels     

                        '3426','3433','3435','3441', # Lag between event and sound is extremely variable and usually >>500ms (something wrong)
                        


                                                    ##### ONES WHOSE FOLDERS REMAIN IN THE MEG DATA FOLDER

                        '2632',  '9000',  #AD          # Too many noisy channels,    
                        '9001', #CH                    # A lot of noise due to people in the room before

                        '9007','9014','9015','9016',#AD # Not recorded from due to metal present on the day
                        '9017',                         # Not recorded from due to problems in the Lab
                        
                        ]






def set_up_participants(num_adults=None,num_children=None):
    ''' 
    Find the participant IDs who will be analysed
    '''
    os.chdir(CHILD_BASE)
    child_participant_strings = []
    base_contents = os.listdir('.')
    for item in base_contents:
        if os.path.isdir(item):
            if len(item) == 4: # All participants have a four-number code
                try:
                    participant_number = int(item)
                    child_participant_strings.append(item)
                except:
                    print("Not a participant: ", item )
    child_participant_strings = sorted(list(set(child_participant_strings)-set(PROBLEM_PARTICIPANTS)))
    child_participant_strings = child_participant_strings[0:num_children] if num_children!=None else child_participant_strings
    #print("Children: ", child_participant_strings)

    # Set up adult participants
    os.chdir(ADULT_BASE)
    adult_participant_strings = []
    base_contents = os.listdir('.')
    for item in base_contents:
        if os.path.isdir(item):
            if len(item) == 4: # All participants have a four-number code
                try:
                    participant_number = int(item)
                    adult_participant_strings.append(item)
                except:
                    print("Not a participant: ", item )

    adult_participant_strings = sorted(list(set(adult_participant_strings)-set(PROBLEM_PARTICIPANTS)))
    adult_participant_strings = adult_participant_strings[0:num_adults] if num_adults!=None else adult_participant_strings
    
    # Enforce that these are lists
    if type(child_participant_strings)!=type([]): # It's a singleton
        child_participant_strings = [child_participant_strings]
    if type(adult_participant_strings)!=type([]): # It's a singleton
        adult_participant_strings = [adult_participant_strings]
    
    return child_participant_strings, adult_participant_strings


child_participant_strings, adult_participant_strings = set_up_participants() # > set_up_participants found in experiment.py
