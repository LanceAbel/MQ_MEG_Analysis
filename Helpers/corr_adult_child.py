import mne
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
from scipy import stats
from scipy.stats.stats import pearsonr
import numpy as np

delay = 3100
child_noisy = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\9013\\9013_HC_22_02_03_EMPTY_ROOM_CHILD_COMPAREvAdult_3300msDelayed.con'
child_filt = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\9013\\9013_HC_22_02_03_EMPTY_ROOM_CHILD_COMPAREvAdult_3300msDelayed_tspca.con'
adult_noisy = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\9013\\9013_HC_22_02_03_EMPTY_ROOM_ADULT_COMPAREvCHILD.con'
adult_filt = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\9013\\9013_HC_22_02_03_EMPTY_ROOM_ADULT_COMPAREvCHILD_tspca.con'

# delay = 2400
# child_noisy = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\9017\\9017_SB_22_02_28_EMPTY_ROOM_CHILD_COMPAREvAdult_2400msDelayed.con'
# child_filt = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\9017\\9017_SB_22_02_28_EMPTY_ROOM_CHILD_COMPAREvAdult_2400msDelayed_tspca.con'
# adult_noisy = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\9017\\9017_SB_22_02_28_EMPTY_ROOM_ADULT_COMPAREvCHILD.con'
# adult_filt = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\9017\\9017_SB_22_02_28_EMPTY_ROOM_ADULT_COMPAREvCHILD_tspca.con'

def get_raw(path_con, trigger_channels):

    raw = mne.io.read_raw_kit(path_con,  # change depending on file i want
    mrk=None,
    elp=None,
    hsp=None,
    stim = trigger_channels,
    slope="+",
    stim_code="channel",
    stimthresh=2, # if self.is_adult_folder else 1,  # 2 for adults
    preload=True,
    verbose=False                    
    )

    return raw


child_noisy = get_raw(path_con=child_noisy,trigger_channels = list(range(145, 145+7)))
#child_filt = get_raw(path_con=child_filt,trigger_channels = list(range(145, 145+7)))
adult_noisy = get_raw(path_con=adult_noisy,trigger_channels = list(range(193, 193+7)))
#adult_filt = get_raw(path_con=adult_filt,trigger_channels = list(range(193, 193+7)))

child_noisy_data = child_noisy.to_data_frame()['MEG 004'].values
#child_filt_data = child_filt.to_data_frame()['MEG 004'].values
adult_noisy_data = adult_noisy.to_data_frame()['MEG 007'].values
#adult_filt_data = adult_filt.to_data_frame()['MEG 007'].values

print("Corr unfiltered without adjusting times", pearsonr(child_noisy_data, adult_noisy_data)[0])
#print("Corr filtered without adjusting times", pearsonr(child_filt_data, adult_filt_data)[0])

plt.scatter(child_noisy_data,adult_noisy_data)
plt.show()
#plt.scatter(child_filt_data,adult_filt_data)
#plt.show()


for time_delay in list(np.arange(delay-2100, delay+2100, 700)):
	adult_noisy_data_delayed = adult_noisy_data[time_delay:]
	#adult_filt_data_delayed = adult_filt_data[time_delay:]
	child_noisy_data_shortened = child_noisy_data[0:-time_delay]
	#child_filt_data_shortened = child_filt_data[0:-time_delay]

	#print(len(child_noisy_data),len(adult_noisy_data),len(adult_noisy_data_delayed),len(child_noisy_data_shortened))
	
	corr_dirty = pearsonr(child_noisy_data_shortened, adult_noisy_data_delayed)[0]
	#corr_filt = pearsonr(child_filt_data_shortened, adult_filt_data_delayed)[0]
	print("Time delay ms %s correlation dirty %s"%(str(time_delay),pearsonr(child_noisy_data_shortened, adult_noisy_data_delayed)[0]))
	#print("Time delay ms %s correlation filt %s"%(str(time_delay),pearsonr(child_filt_data_shortened, adult_filt_data_delayed)[0]))
	
	#plt.scatter(child_noisy_data_shortened,adult_noisy_data_delayed)
	x = child_noisy_data_shortened
	y = adult_noisy_data_delayed
	

	#regression part
	# method A
	# m, b = np.polyfit(x, y, 1) # m = slope, b=intercept.
	# plt.plot(x,y,'.')
	# plt.plot(x, m*x + b)    #add line of best fit.
	# method B
	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	line = slope*x+intercept
	plt.plot(x, line, 'r', label='y={:.2f}x+{:.2f}'.format(slope,intercept))
	#plt.plot(x, line, 'r', label='corr={:.2f}'.format(corr_dirty))
	plt.scatter(x,y, color="k", s=3.5)
	plt.legend(fontsize=9)
	plt.title("Correlation at delay %s ms"%(str(time_delay-delay)))
	plt.xlabel("Child sensor data")
	plt.ylabel("Adult sensor data")
	plt.show()
	#plt.scatter(child_filt_data_shortened,adult_filt_data_delayed)
	#plt.show()