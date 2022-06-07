## Other helper functions
# Find files in a folder
from stat import *
from config import *
#from participant_data import *
# Importing the library
import psutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure




#OS   

def cpu_used(time_interval_secs=None):
    '''Returns fraction of CPU used over a timeframe'''
    return psutil.cpu_percent(time_interval_secs)

print("Cpu Used as %", cpu_used(time_interval_secs=2))


def free_memory():
    total_system_memory_excl_swap = psutil.virtual_memory()[0]/(2**30)
    percent_memory_used = psutil.virtual_memory()[2]
    free_memory_gb = psutil.virtual_memory()[4]/(2**30)
    return free_memory_gb
print("Memory Free GB", free_memory())


def wait_until_memory_free(required_memory = 6, max_wait_time_mins = 5):
    '''Helps pause a thread until enough memory is free
    Will re-attempt if no memory freed up over max_wait_time_mins
    
    Designed to allow us to run multiple threads whilst preventing a MemoryError
    ''' 

    memory_free = False
    secs_paused_for = 0
    pause_time = 10 # secs
    while not memory_free:
        mem_free = free_memory()
        if mem_free > required_memory:
            memory_free = True
            show_paused = False
        else:
            if secs_paused_for/60 < max_wait_time_mins: # End after waiting 5 mins unsuccessfully for memory to free up
                print("Pausing for 10 secs, memory free only ", mem_free, " GB. Mins paused for ", secs_paused_for/60)
                secs_paused_for+=pause_time
                time.sleep(pause_time)
            else:
                print("Gave up waiting for memory to become free after ", secs_paused_for/60, " mins ")
                memory_free = True


def walktree(top, callback):
    '''recursively descend the directory tree rooted at top,
       calling the callback function for each regular file'''

    for f in os.listdir(top):
        pathname = os.path.join(top, f)
        mode = os.stat(pathname)[ST_MODE]
        if S_ISDIR(mode):
            # It's a directory, recurse into it
            walktree(pathname, callback)
        elif S_ISREG(mode):
            # It's a file, call the callback function
            callback(pathname)
        else:
            # Unknown file type, print a message
            #print('Skipping %s' % pathname)
            pass 


def visitfile(file):
    global all_files
    all_files.append(file)
    
def find_all_files(folder):
    global all_files
    all_files = []
    walktree(folder, visitfile)
    return all_files



# Logging 
num_calls = {}
def logger(stack, extras=None):
    '''Used to track which function is being run and where Exceptions arise, if any'''
    if LOGGING:
        try:
            stack = stack[0][3]
            if stack in num_calls:
                num_calls[stack] = num_calls[stack] +1
            else:
                num_calls[stack] = 1
            printstr = "####### LOGGER:  Stack %s , call # %s for this participant"%(str(stack),str(num_calls[stack]))
            if extras!=None:
                printstr += " Extra info %s"%(str(extras))
            print(printstr)
        except:
            print("PROBLEM LOGGING, ", len(stack))
            print("PROBLEM LOGGING, ", len(stack[0]))
            

# Data manipulation
def try_add(dct,key,value):
    
    if key in dct.keys():
        #print(key,dct[key])
        dct[key] = dct[key] + [value]
    else:
        dct[key] = [value]   


def try_return(obj,attr):
    if hasattr(obj,attr):
        return getattr(obj,attr)
    else:
        return None
           
def try_del(items):
    '''Try to delete a variable'''
    for item in items:
        if item in globals() or item in locals():
            print("Deleting ", item)
            try:
                del item
            except:
                pass    

def try_del_attr(obj, attr):
    '''Try to delete an attribute of an object'''    
    if hasattr(obj,attr):
        delattr(obj, attr)
        
def del_attributes(pctp,attrs_to_delete):
    # Delete information in each participant that is not necessary to do group analysis
    for attrib in attrs_to_delete:
        try_del_attr(pctp, attrib)
    gc.collect()   
 

#Data manip
def split_list(a_list):
    half = len(a_list)//2
    return a_list[:half], a_list[half:]

    
def equalObs(x, nbin):
    '''Make nbins of equal size from data x'''
    nlen = len(x)
    return np.interp(np.linspace(0, nlen, nbin + 1),
                     np.arange(nlen),
                     np.sort(x))


def absolute(x):
    return abs(x)
    
def find_cutoff(dct,percentage_cutoff=80):
    '''
    Finds the smallest key value in a dictionary of key values for which the sum of the dictionary's values for keys <= itself is greater than the percentage_cutoff
    '''
    keys = sorted(dct.keys())
    sum_vals = sum(dct.values())
    cumulative_total = 0
    cutoff_found = False
    for key in keys:
        if cutoff_found:
            cutoff_key = key
            return cutoff_key
        else:
            cumulative_total+=dct[key]
            if cumulative_total/sum_vals > percentage_cutoff/100:
                cutoff_found = True
                    
def get_max_keys(dct):
    # Get two most frequently identified ICA indices to exclude
    max_key = max(dct, key= dct.get) if len( dct) > 0 else None
    second_largest_key = None
    ica_counts = dct.values()
    if len(dct) > 1:
        second_largest = sorted(ica_counts)[max(-2,-len( dct))]
    else:
        second_largest = None

    for count_num in ica_counts:
        for ica_index in dct.keys():
            if dct[ica_index] == second_largest:
                second_largest_key = ica_index
    #print("Max ", max_key, " 2nd: ", second_largest_key)
    return(max_key,second_largest_key)

def increment_ct(dct,key):
    if key in dct.keys():
        dct[key] = dct[key] + 1
    else:
        dct[key] = 1

 
def RMS_DF(df):
    '''Return the root mean square of all the values across MEG data columns in a df for each time period'''
    times = df['time'].values
    indices = df.index.values
    rms_values = []
    for index in indices: # For each time
        sum_sq = 0
        num_cols = 0
        for col in df.columns.values:
            if "MEG" in col:
                sum_sq+=df[col][index]**2
                num_cols+=1
        mean_sq = sum_sq/num_cols
        rms = math.sqrt(mean_sq)
        rms_values.append(rms)
    return rms_values

        
def sum_df(df,mode='plain'):
    df['sum'] = [0 for x in range(0,len(df.index.values))]
    for col in df:
        if "MEG" in col and "time" not in col:
            if mode=="plain":
                df['sum'] = df['sum'] + df[col]       
            elif mode == "abs":
                df[col+"abs"] = abs(df[col])
                df['sum'] = df['sum'] + df[col+'abs']
                df = df.drop(col+'abs', 1)
            elif mode == "square":
                df[col+"sq"] = df[col]**2
                df['sum'] = df['sum'] + df[col+'sq']
                df = df.drop(col+'sq', 1)
    return df

def avg_df(df,mode='plain'):
    df = sum_df(mode=mode)
    cols = [chan for chan in df.columns.values if "MEG" in chan]
    L = len(cols)
    df['avg'] = df['sum']*1.0/L
    return df


# Alter MEG sensor info
def centre_sensor_locations(evoked):
    evoked_copy = copy.deepcopy(evoked)
    ch_names = evoked_copy.info['ch_names']
    C  = 1
    x_s = []
    y_s = []
    z_s = []
    others = {}
    for ch in ch_names:
        for dct in evoked_copy.info['chs']:
            if dct['ch_name'] == ch:
                loc = dct['loc']
                x_s.append(loc[0])
                y_s.append(loc[1])
                z_s.append(loc[2])
#                 for r in range(3,11): # Other dimensions don't do anything
#                     if r in others.keys():
#                         others[r].append(loc[r])
#                     else:
#                         others[r] = [loc[r]]

    avg_ = {}
    avg_x = sum(x_s)/len(x_s)
    avg_y = sum(y_s)/len(y_s)
    avg_z = sum(z_s)/len(z_s)
    for r in range(3,11): # Other dimensions don't do anything
        avg_[r] = sum(others[r])/len(others[r])

    #sq_before = (avg_x**2+avg_y**2+avg_z**2)**0.5
    for ch in ch_names:
        r = 0
        for dct in evoked.info['chs']:
            if dct['ch_name'] == ch:
                evoked_copy.info['chs'][r]['loc'][0] = (evoked_copy.info['chs'][r]['loc'][0]-avg_x)/C
                evoked_copy.info['chs'][r]['loc'][1] = (evoked_copy.info['chs'][r]['loc'][1]-avg_y)/C
                # Z dimension screws things up (polar co-ordinates?)
                # evoked_copy.info['chs'][r]['loc'][2] = (evoked_copy.info['chs'][r]['loc'][2]-avg_z)/C 
                # for r in range(3,11): # Other dimensions don't do anything
                #     evoked_copy.info['chs'][r]['loc'][r] = (evoked_copy.info['chs'][r]['loc'][r]+avg_[r])/C   
            r+=1

    return evoked_copy, avg_x, avg_y, avg_z


def find_radius(evoked, radius_=None):
    if radius_:
        scalar=radius_*1.07/0.26824116413966   #26.8cm child helmet 
    else:
        scalar=radius*1.07/0.26824116413966
    
    r = 0
    max_dist = -9999
    for ch in evoked.info['chs']:
        loc = ch['loc']
        x,y,z = loc[0],loc[1], loc[2]
        for ch_inner in evoked.info['chs']:
            loc_inner = ch_inner['loc']
            x_inner,y_inner,z_inner = loc_inner[0],loc_inner[1], loc_inner[2]    
            dist = ((x_inner-x)**2 + (y_inner-y)**2 + (z_inner-z)**2)**0.5
            if dist > max_dist:
                max_dist = dist
        r+=1
    return scalar*max_dist


def align_channels_paul_individual(ptcp):

    data_dir = "E:\Downloads\Change_channel_mapping\\"
    MEG125_raw = mne.io.read_raw_kit(data_dir + "MEG125_sample.con",
                                        allow_unknown_format=False,
                                        verbose=True)

    data_dir = "E:\Downloads\Change_channel_mapping\\"
    # load up the mapping lists - stored as .mat files
    MEG125 = scipy.io.loadmat(
        data_dir + "labels_systemmatch_child_sensors.mat", simplify_cells=True)
    MEG160 = scipy.io.loadmat(
        data_dir + "labels_systemmatch_adult_sensors.mat", simplify_cells=True)    
    MEG125_picks = MEG125["lab125"]
    MEG160_picks = MEG160["lab160"]

    # Change channel names
    for idx, ch in enumerate(MEG160_picks):
        MEG160_picks[idx] = ch.replace("AG", "MEG ")
    for idx, ch in enumerate(MEG125_picks):
        MEG125_picks[idx] = ch.replace("AG", "MEG ") 


    MEG125_picks = list(set(MEG125_picks)-set(['MEG 072'])) # 072 is missing from mapping
    MEG160_picks = list(set(MEG160_picks)-set(['MEG 072']))    
    for m in range(126,161):
        MEG160_picks = list(set(MEG160_picks)-set(["MEG "+str(m).zfill(3)]))
    MEG125_raw.pick_channels(MEG125_picks)    


    if ptcp.is_adult_system:
        ptcp.get_raw(trigger_channels_known=False)
        events = mne.find_events(
                                ptcp.raw,
                                output="onset",
                                consecutive=False,
                                min_duration=0,
                                shortest_event=1,  # 5 for adults
                                mask=None,
                                uint_cast=False,
                                mask_type="and",
                                initial_event=False,
                                verbose=None
                                )

        #MEG160_epochs = mne.Epochs(MEG160_raw, events, tmin=-0.1, tmax=0.4, preload=True)
        #MEG160_epochs.pick_types("mag")

        #info = mne.create_info(ch_names, self.epochs_ransac_autoreject.info['sfreq'], ch_types='mag', verbose=True)
        for attrib in ['raw_cleaned','evoked_generic','evoked_all','epochs','epochs_ransac','epochs_ransac_autoreject','epochs_ransac_autoreject_unstandardised']:
            if hasattr(ptcp, attrib):
                obj = getattr(ptcp, attrib)
                if type(obj) == type({}):
                    keys = list(obj.keys())
                    for key in keys:
                        try:
                            obj[key].pick_types("mag")
                        except:
                            print("PROBLEM PICKING MAG")
                            print(obj[key].info['ch_names'])
                        try:
                            obj[key].pick_channels(MEG160_picks)
                        except:
                            print("PROBLEM PICKING CHANNELS")
                            print(MEG160_picks)
                            for ch in MEG160_picks:
                                if ch not in obj[key].info['ch_names']:
                                    print("Mismatch 1", ch)
                            for ch in obj[key].info['ch_names']:                                    
                                if ch not in MEG160_picks:
                                    print("Mismatch 2", ch)
                        obj[key].reorder_channels(MEG160_picks)
                        obj[key].info["chs"] = MEG125_raw.info["chs"]
                        obj[key].info["ch_names"] = MEG125_raw.info["ch_names"]
                else:
                    obj.pick_types("mag")
                    obj.pick_channels(MEG160_picks)
                    obj.reorder_channels(MEG160_picks)
                    obj.info["chs"] = MEG125_raw.info["chs"]
                    obj.info["ch_names"] = MEG125_raw.info["ch_names"]
        delattr(ptcp, 'raw')



# Stats
def fisher_r_to_z(r):
    # fisher r to z transformation
    return 0.5*math.log((1+r)/(1-r)) # fisher r to z transformation
def fisher_z_to_r(z):
    # fisher z to r transformation
    return ( (math.e)**(2*z)-1)/(  (math.e)**(2*z)+1)

def compare_two_corrs(r1,r2,n1,n2,critical_z = 1.96):
    '''
    Returns True if the correlations are significantly different (at a p=0.05 level), otherwise False
    '''
    
    z = (fisher_r_to_z(r1) -  fisher_r_to_z(r2))/math.sqrt((1/(n1-3)) + (1/(n2-3)))
    if z < -critical_z or z> critical_z:
        return True
    else:
        return False

def corr_coeff_t_value(n,corr):
    # Test the significance of the correlation co-efficient
    t = corr/math.sqrt((1-corr**2)/(n-2))
    return t

def calc_z_score(col):
    #assert(np.std(col)!=0),"Cannot standardise, standard deviation is zero "+str(col)
    if np.std(col)==0:
        print("WARN: Std of column is zero")
        return col
    else:
        #print(np.std(col))
        return (col-np.mean(col))/np.std(col) /(10**13) # Divider to keep the signal roughly the same in value

def standardise(col):
    #assert(np.std(col)!=0),"Cannot standardise, standard deviation is zero "+str(col)
    if np.std(col)==0:
        print("WARN: Std of column is zero")        
        return col
    else:    
        diff = (col-np.mean(col))*1.0
        return diff/np.std(col)

def find_scalar(df_1,df_2):
    '''After taking a z score of the signals in df_2, multiply all values by some number to give a signal of approx the same magnitude'''
    scalars=[]
    for col in df_1.columns.values:
        if "MEG" in col:
            scalar = np.std(df_1[col])/np.std(df_2[col])
            scalars.append(scalar)
    multiplier = sum(scalars)/len(scalars)
    print(multiplier)
    return multiplier   

def apply_z_score_new(ptcp, STANDARDISE=True, channel_wise=False):
    '''This is placed here in case there is an old Participant object without this method'''
    if STANDARDISE:
        for attrib in ['raw_cleaned','evoked_generic','evoked_all','epochs','epochs_ransac','epochs_ransac_autoreject']:
            if hasattr(ptcp,attrib):
                obj = getattr(ptcp, attrib)
                if channel_wise == True:
                    if type(obj) == type({}):
                        keys = list(obj.keys())
                        for key in keys:
                            print("Applying 1, ", key)
                            obj[key].apply_function(fun=calc_z_score, n_jobs=4)
                    else:
                        obj.apply_function(fun=calc_z_score, n_jobs=4)
                else:
                    try:
                        obj.apply_function(fun=calc_z_score, n_jobs=4)#, channel_wise=False)  
                    except Exception as e:
                        keys = list(obj.keys())
                        for key in keys:
                            print("Applying 2, ", key)
        ptcp.signals_of_interest()
        ptcp.applied_z_score = True
    else:
        ptcp.applied_z_score = False


def colour_regression(x,y,title,xlabel,ylabel,order=1):
    plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True
    #fig, ax = plt.subplots()
    _ = plt.scatter(x, y, c=x, cmap='plasma')
    z = np.polyfit(x, y, order)
    p = np.poly1d(z)
    plt.plot(x, p(x), "r-o")
    print(title,type(title))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
        
def independent_ttest(data1, data2, alpha):
    # calculate means
    mean1, mean2 = np.mean(data1), np.mean(data2)
    # calculate standard errors
    se1, se2 = sem(data1),  sem(data2)
    # standard error on the difference between the samples
    sed = math.sqrt(se1**2.0 + se2**2.0)
    # calculate the t statistic
    t_stat = (mean1 - mean2) / sed
    # degrees of freedom
    df = len(data1) + len(data2) - 2
    # calculate the critical value
    cv = t.ppf(1.0 - alpha, df)
    # calculate the p-value
    p = (1.0 - t.cdf(abs(t_stat), df)) * 2.0
    # return everything
    return t_stat, df, cv, p

# function for calculating the t-test for two dependent samples
def dependent_ttest(data1, data2, alpha):
    # calculate means
    mean1, mean2 = mean(data1), mean(data2)
    # number of paired samples
    n = len(data1)
    # sum squared difference between observations
    d1 = sum([(data1[i]-data2[i])**2 for i in range(n)])
    # sum difference between observations
    d2 = sum([data1[i]-data2[i] for i in range(n)])
    # standard deviation of the difference between means
    sd = sqrt((d1 - (d2**2 / n)) / (n - 1))
    # standard error of the difference between the means
    sed = sd / sqrt(n)
    # calculate the t statistic
    t_stat = (mean1 - mean2) / sed
    # degrees of freedom
    df = n - 1
    # calculate the critical value
    cv = t.ppf(1.0 - alpha, df)
    # calculate the p-value
    p = (1.0 - t.cdf(abs(t_stat), df)) * 2.0
    # return everything
    return t_stat, df, cv, p


##### Functions that  plot the correlation over different channels,   X time
def predictor_stat_chan_time(stats,title):
    '''
    Plot a colour map of a statisic (e.g. correlation, between HGF and a MMR) across time and channel space
    Input is a one-person statistic X time dictionary
    '''
    
    times_ref = stats.keys()
    
    a = []
    for time in times_ref:
        a.append(list(stats[time].values()))


    ticks = list(stats.keys())
    fig = plt.figure(figsize=(9, 4.2))
    plt.figure(figsize=(9, 4.2))
    plt.res=64*10
    plt.tcbar_fmt='%5.3f'
    
    # We need to draw the canvas, otherwise the labels won't be positioned and 
    # won't have values yet.
    fig.canvas.draw()

    plt.imshow(a, cmap='hot', interpolation='nearest')
    
    plt.ylabel("ms in EPOCH")
    plt.xlabel("Channel #")
    gap = max(times_ref) - min(times_ref) # EPOCH_END_ANALYSIS-EPOCH_START_ANALYSIS
    labels_y = [-60+min(times_ref) + (int(gap)/(2-1)) * item for item in range(0,5)]
    #labels_y = [min(times_ref) + (int(gap)/(5-1)) * item for item in range(0,5)]
   
    print(labels_y)
    ax = fig.add_subplot(111)  
    ax.set_title(title)      
    ax.set_yticks(labels_y)
    ax.set_yticklabels(labels_y)
    plt.title([title])
    plt.show()


def predictor_stat_chan_time_new(stats_all_ptcps,statistic,ptcp_string,title):
    '''
    Plot a colour map of a statistic (e.g. correlation) between a predictor and a MMR, across time and channel space
    Input is a multi-person statistic X participant X time dictionary
    '''
    times_ref = stats_all_ptcps[statistic][ptcp_string].keys()
    
    a = []
    for time in times_ref:
        a.append(list(stats_all_ptcps[statistic][ptcp_string][time].values()))


    fig = plt.figure(figsize=(9, 4.2))
    plt.figure(figsize=(9, 4.2))
    plt.res=64*10
    #plt.size=2
    plt.tcbar_fmt='%5.3f'
    #fig.canvas.draw()    
    plt.imshow(a, cmap='hot', interpolation='nearest')
    
    
    ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.set_aspect('equal')

    cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    plt.colorbar(orientation='horizontal',label=statistic)
    #plt.clim(-1, 1)
    
    plt.ylabel("ms in EPOCH")
    plt.xlabel("Channel #")
    gap = max(times_ref) - min(times_ref) # EPOCH_END_ANALYSIS-EPOCH_START_ANALYSIS
    labels_y = [-60+min(times_ref) + (int(gap)/(2-1)) * item for item in range(0,5)]
    #labels_y = [min(times_ref) + (int(gap)/(5-1)) * item for item in range(0,5)]
    #print(labels_y)
    ax.set_yticklabels(labels_y)
#                                 vmin=0,cmap='RdBu_r',
#                             res=64*3, size=2, cbar_fmt='%5.3f',
            
    plt.show()






 