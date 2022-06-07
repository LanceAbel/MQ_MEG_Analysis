
windows_machine = ispc 
if windows_machine
    % FT toolbox
    addpath(genpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Paul, Jordan\Jordan\Jordan Code\"))
    addpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Matlab_toolboxes\fieldtrip-master\fieldtrip-master")
    ft_defaults
    addpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Matlab_toolboxes\fieldtrip-master\fieldtrip-master\external\eeglab")
    addpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Matlab_toolboxes\fieldtrip-master\fieldtrip-master\external\mffmatlabio\private")
    %SBL
    addpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\SequentialBayesianLearning-master\Sequential_Bayesian_learner-master")
    %SPM
    addpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\spm12-master\spm12-master")
    % EEG lab
    addpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Matlab_toolboxes\eeglab\functions\popfunc")
    addpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Matlab_toolboxes\eeglab\functions\adminfunc")
    addpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Matlab_toolboxes\eeglab\functions\sigprocfunc")
    addpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Matlab_toolboxes\eeglab\plugins\firfilt2.4")
    % Roving MMN
    addpath(genpath("E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Roving_MMN_2020-master\Roving_MMN_2020-master\preprocessing_scripts_ME125_phase2_YSun\"))
    addpath(genpath("E:\oDrive\OneDrive\Docs Sync\MyProjects\Code\Matlab Code\"))
    
    %cd("E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Child_MEG\2872\")
else
    % FT toolbox
    addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Matlab_toolboxes/fieldtrip-master/fieldtrip-master')
    ft_defaults
    addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Matlab_toolboxes/fieldtrip-master/fieldtrip-master/external/eeglab')
    addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Matlab_toolboxes/fieldtrip-master/fieldtrip-master/external/mffmatlabio/private')
    %SBL
    addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/SequentialBayesianLearning-master/Sequential_Bayesian_learner-master')
    %SPM
    addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/spm12-master/spm12-master')
    % EEG lab
    addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Matlab_toolboxes/eeglab/functions/popfunc')
    addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Matlab_toolboxes/eeglab/functions/adminfunc')
    addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Matlab_toolboxes/eeglab/functions/sigprocfunc')
    addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Matlab_toolboxes/eeglab/plugins/firfilt2.4')
    % Roving MMN
    addpath(genpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Roving_MMN_2020-master/Roving_MMN_2020-master/preprocessing_scripts_ME125_phase2_YSun'))
    % addpath(genpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Roving_MMN_2020-master/Roving_MMN_2020-master/preprocessing_scripts_ME125_phase2_YSun/Others'))
    % Or sub-paths - had to add for MAC run
    %addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Roving_MMN_2020-master/Roving_MMN_2020-master/preprocessing_scripts_ME125_phase2_YSun/Others')
    %addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Roving_MMN_2020-master/Roving_MMN_2020-master/preprocessing_scripts_ME125_phase2_YSun/Group Level')
    %addpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Roving_MMN_2020-master/Roving_MMN_2020-master/preprocessing_scripts_ME125_phase2_YSun/Group Level/Stats')    
    
    addpath(genpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Macquarie-MEG-Research/MQ_MEG_Scripts-master/MQ_MEG_Scripts-master/Tools'))
    addpath(genpath('/Users/user/odrive/OneDrive/Docs Sync/Jobs And Money/Careers/Study/Uni/MAC/MQ MRes/Coding/Macquarie-MEG-Research/MQ_MEG_Scripts-master/MQ_MEG_Scripts-master/Preprocessing'))
    %cd('/Users/user/Downloads/MEG Temp/2872')
    
end

