# Deletes previously-generated results files
OLDFILE_CUTOFF_DAYS = 0.5 # Delete files older than this (in days)


import os
import glob

timestamp_file = r'C:\Users\Lance\Desktop\CanDelete\\Temp.txt'
with open(timestamp_file, 'w') as f:
    f.write('Temp file')
time_now = os.path.getmtime(timestamp_file)
EXP_BASE          = r'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Experiments\\'
os.chdir(EXP_BASE)
base_contents = os.listdir('.')
print(base_contents)
num_files = 0
for x in base_contents[2:]:
    subdir = EXP_BASE+x+"\\"
    if os.path.isdir(subdir):
        os.chdir(subdir)
        print(subdir)
        csvs = glob.glob(f"{subdir}*.csv")        
        txts = glob.glob(f"{subdir}*.txt")
        jpgs = glob.glob(f"{subdir}*.jpg")
        jpegs = glob.glob(f"{subdir}*.jpeg")
        for file in csvs+txts+jpgs+jpegs:
            time_mod = os.path.getmtime(file)
            days_since = (time_now - time_mod)/60/60/24
            if days_since > OLDFILE_CUTOFF_DAYS: # File wasn't found or is old
                print("Deleting ", file, days_since*24*60, "mins")
                os.remove(file)
                num_files+=1

    subdir = EXP_BASE+x+"\\Standardised\\"
    if os.path.isdir(subdir):
        os.chdir(subdir)
        print(subdir)
        csvs = glob.glob(f"{subdir}*.csv")        
        txts = glob.glob(f"{subdir}*.txt")
        jpgs = glob.glob(f"{subdir}*.jpg")
        jpegs = glob.glob(f"{subdir}*.jpeg")
        for file in csvs+txts+jpgs+jpegs:
            time_mod = os.path.getmtime(file)
            days_since = (time_now - time_mod)/60/60/24
            if days_since > OLDFILE_CUTOFF_DAYS: # File wasn't found or is old
                print("Deleting ", file, days_since*24*60, "mins")
                os.remove(file)
                num_files+=1

print("Deleted ", num_files, " files")