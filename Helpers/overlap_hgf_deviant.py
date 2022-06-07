## Calculate correlations between predictors

import matplotlib.pyplot as plt
import numpy as np
import pprint

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
def list_at_indices(list, indices):
	return [list[index] for index in indices if index<=len(list)]

tones_file	= r"E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\\min 1, max 7, skewed reps uniform tones.txt"
tones = get_data(tones_file)
tones = process_file(tones)
deviants = [0] + [(tones[x]!=tones[x-1])*1.0 for x in range(1,len(tones))]

indices_deviants = [x for x in range(0,len(deviants)) if deviants[x] == 1]
print(indices_deviants)




## CORRELATION TO DEVIANT FILE
# HGFs
vanilla_folder = r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\HGF Original\\'
#mod_folder = r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\HGF Settings999\\'
mod_folder = r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\HGF Settings Used\\'
applied_folder = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Common\\'
for hgf_size in ['1730']:
#for hgf_size in [int(346*x) for x in range(1,6)] + [2000]:
	# Vanilla HGF
	#pred_file = vanilla_folder+'HGF'+str(hgf_size)+' PE2_new.txt'	

	# Modified HGF
	#pred_file = mod_folder+'HGF'+str(hgf_size)+' PE2_mod_new.txt'
	pred_file =  mod_folder+ 'HGF'+str(hgf_size)+' PE2_mod_baked_new.txt'

	hgf = get_data(pred_file)
	hgf = process_file(hgf)

	size_to_analyse = min(int(hgf_size),len(hgf),len(deviants))
	print(size_to_analyse)
	print(np.corrcoef(deviants[0:size_to_analyse],hgf[0:size_to_analyse] ) [0] [1])



# # Corr of Other predictors to deviant/predeviant
# deviants_mod = deviants
# for pred_file in [ applied_folder+'BS1730_new.txt',
# 					applied_folder+'PS1730_new.txt',
# 					applied_folder+'CS1730_new.txt',
# 					#r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\Other Surprise\Raw Probability\\min 1, max 7, skewed reps uniform tones_identities PROB.txt'
# 					]:
# 	#print(pred_file)
# 	pred = get_data(pred_file)
# 	try:
# 		pred = process_file(pred)	
# 	except:
# 		pred = process_file2(pred)
# 	size_to_analyse = min(len(pred),len(deviants_mod))
# 	print(size_to_analyse)
# 	for r in range(0,30): # Show the actual values
# 		pass
# 		#print(r, deviants_mod[r], pred[r])
# 	print(pred_file, np.corrcoef(deviants_mod[0:size_to_analyse],pred[0:size_to_analyse] ) [0] [1])					




## CORRELATION TO EACH OTHER ON DEVIANTS
deviants_seq_at_deviants = list_at_indices(deviants,indices_deviants)	
corr_matrix = {}
m = 0
print("@@@@@@@@@@@@")
print([int(346*x) for x in range(1,6)] + [2000])
#for hgf_size in ['1730']:
for hgf_size in [int(346*x) for x in range(1,6)]:# + [2000]:
	print("@@@@@@@ HGF SIZE @@@@@@@@ " +str(hgf_size))
	pred_files = 	{
					"VanillaPE2": vanilla_folder+'HGF'+str(hgf_size)+' PE2_new.txt',
					 "VanillaPWPE2": vanilla_folder+'HGF'+str(hgf_size)+' PWPE2_new.txt',
					 "VanillaPE3": vanilla_folder+'HGF'+str(hgf_size)+' PE3_new.txt',
					 "VanillaPWPE3": vanilla_folder+'HGF'+str(hgf_size)+' PWPE3_new.txt',
					# "ModifiedPE2": mod_folder+'HGF'+str(hgf_size)+' PE2_mod_new.txt',
					# "ModifiedPWPE2": mod_folder+'HGF'+str(hgf_size)+' PWPE2_mod_new.txt',	
					# "ModifiedPE3": mod_folder+'HGF'+str(hgf_size)+' PE3_mod_new.txt',
					# "ModifiedPWPE3": mod_folder+'HGF'+str(hgf_size)+' PWPE3_mod_new.txt',														
					"ModifiedBakedPE2": mod_folder+'HGF'+str(hgf_size)+' PE2_mod_baked v2_new.txt',
					"ModifiedBakedPE3": mod_folder+'HGF'+str(hgf_size)+' PE3_mod_baked v2_new.txt',		
					"ModifiedBakedPWPE2": mod_folder+'HGF'+str(hgf_size)+' PWPE2_mod_baked v2_new.txt',
					"ModifiedBakedPWPE3": mod_folder+'HGF'+str(hgf_size)+' PWPE3_mod_baked v2_new.txt'										
					}

	for pred in pred_files.keys():
		for pred_ in pred_files.keys():

			if pred != pred_:

				if hgf_size not in corr_matrix.keys():
					corr_matrix[hgf_size] = {}
				

				hgf_1 = get_data(pred_files[pred])
				hgf_1 = process_file(hgf_1)
				hgf_2 = get_data(pred_files[pred_])
				hgf_2 = process_file(hgf_2)
				size_to_analyse = min(len(hgf_1),len(hgf_2))				
				# if m < 10:
				# 	print(m, hgf_size, pred)
				# 	plt.plot(hgf_1[0:size_to_analyse],linestyle="",marker="o")
				# 	plt.title(pred+ " Length: "+str(size_to_analyse)) # hgf_size
				# 	plt.legend([pred,size_to_analyse])
				# 	plt.show()
				m+=1

				

				# hgf_1 = list_at_indices(hgf_1,indices_deviants)
				# hgf_2 = list_at_indices(hgf_2,indices_deviants)				


				size_to_analyse = min(len(hgf_1),len(hgf_2))
				#print("Sizes to analyse ", size_to_analyse)

				corr = np.corrcoef(hgf_1[0:size_to_analyse],hgf_2[0:size_to_analyse] ) [0] [1]
				#print("Correlation b/w ", pred + " and ", pred_, corr)


				if (pred_, pred) not in corr_matrix[hgf_size].keys():
					corr_matrix[hgf_size][pred,pred_] = corr

				#corr2 = np.corrcoef(hgf_1[0:size_to_analyse_dev],deviants_seq_at_deviants[0:size_to_analyse_dev] ) [0] [1]
				#print("Correlation b/w deviant and ", pred, corr2) # will be NAN as the deviants are all == 1
				# if ("deviant", pred) not in corr_matrix[hgf_size].keys():
				# 	corr_matrix[hgf_size][pred,"deviant"] = corr2


for time_hor in corr_matrix.keys():
	print(time_hor)
	for mod1, mod2 in corr_matrix[time_hor].keys():
		print(mod1, mod2, corr_matrix[time_hor][mod1, mod2])

#print(corr_matrix)




## Calculate correlation between HGFs with different settings
hgf_settings_2_predictions = r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\HGF Settings2\\min 1, max 7, skewed reps skewed tones HGF_Settings2_1720 PE2.txt'
hgf_settings_3_predictions = r'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\HGF Settings3\\min 1, max 7, skewed reps skewed tones HGF_Settings3_1720 PE2.txt'
settings_2_preds = get_data(hgf_settings_2_predictions)
settings_2_preds = [float(x.replace("\n","")) for x in settings_2_preds]
settings_3_preds = get_data(hgf_settings_3_predictions)
settings_3_preds = [float(x.replace("\n","")) for x in settings_3_preds]
corr = np.corrcoef(settings_2_preds, settings_3_preds) [0] [1]
print("Correlation between HGFs with different settings ", corr)



# ## Calculate correlation between predictors fitted to data of different length
# for length in corr_matrix.keys():
# 	for length_inner in corr_matrix.keys():
# 		if length!=length_inner:
# 			for key in pred_files.keys():
# 				print(key)
# 				file_1 = pred_files[key]
# 				file_2 = pred_files[key]
# 				file_1 = file_1.replace(str(hgf_size),str(length))
# 				#print(file_1)
# 				file_2 = file_2.replace(str(hgf_size),str(length_inner))
# 				#print(file_2)
# 				hgf_1 = get_data(file_1)
# 				hgf_1 = process_file(hgf_1)
# 				hgf_2 = get_data(file_2)
# 				hgf_2 = process_file(hgf_2)

# 				size_to_analyse = min(len(hgf_1),len(hgf_2))	

# 				corr = np.corrcoef(hgf_1[0:size_to_analyse],hgf_2[0:size_to_analyse] ) [0] [1]
# 				print(length,length_inner,corr)



# ## Calculate correlation between predictor1 and deviant (predictor2) fitted to data of same length
# for length in corr_matrix.keys():
# 	for length_inner in corr_matrix.keys():
# 		if length==length_inner:
# 			for key in pred_files.keys():
# 				if key=="ModifiedBakedPE2":#"VanillaPE2":
# 					print(key)
# 					file_1 = pred_files[key]
# 					file_2 = pred_files[key]
# 					file_1 = file_1.replace(str(hgf_size),str(length))
# 					#print(file_1)
# 					hgf_1 = get_data(file_1)
# 					hgf_1 = process_file(hgf_1)
# 					hgf_1 = hgf_1[1:]
					
# 					hgf_2 = deviants[:-1]

# 					#print(hgf_1[0:30])
# 					#print(hgf_2[0:30])

# 					size_to_analyse = min(len(hgf_1),len(hgf_2))	

# 					corr = np.corrcoef(hgf_1[0:size_to_analyse],hgf_2[0:size_to_analyse] ) [0] [1]
# 					print(length,length_inner,corr)


## Calculate correlation between predictor1 and same predictor fitted to data of different length
for length in corr_matrix.keys():
	for length_inner in corr_matrix.keys():
		if length!=length_inner:
			for key in pred_files.keys():
				if key=="VanillaPE2":
					print(key)
					file_1 = pred_files[key]
					file_2 = pred_files[key]
					file_1 = file_1.replace(str(hgf_size),str(length))
					#print(file_1)
					hgf_1 = get_data(file_1)
					hgf_1 = process_file(hgf_1)
					

					file_2 = pred_files[key]
					file_2 = pred_files[key]
					file_2 = file_2.replace(str(hgf_size),str(length_inner))
					#print(file_1)
					hgf_2 = get_data(file_2)
					hgf_2 = process_file(hgf_2)
	

					#print(hgf_1[0:30])
					#print(hgf_2[0:30])

					size_to_analyse = min(len(hgf_1),len(hgf_2))	

					corr = np.corrcoef(hgf_1[0:size_to_analyse],hgf_2[0:size_to_analyse] ) [0] [1]
					print(length,length_inner,corr)