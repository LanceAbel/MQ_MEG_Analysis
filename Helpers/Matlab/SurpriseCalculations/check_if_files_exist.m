% Checks if any HGF predictor files of length x are missing, outputs them to 
% a) HGFxxxx for adults
% b) HGFxxxx_children for children
% Also checks if any of the Sound Delay.txt files are missing

% Children: 2699 is just a bit shorter.     2724 and 2810, there are no sounds for and were excluded
% Adults: 2607 had 11 repetitions.          2552 and 2689, there are no sounds for and were excluded
NUM_TONES_LONG = 692
%SUBJECTS_TO_RUN = ["", "9002"] % {'9002','9001'} % Comment out line if you want to run on everyone

windows_machine = ispc
if windows_machine
    details_file = 'E:\subject.txt';
else
    details_file = '/subject.txt'
end
fid=fopen(details_file);
tline = fgetl(fid);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);
matlab_general_code = char(tlines(6));
code_folder = char(tlines(7));
code_folder_mac = char(tlines(8));
if ~windows_machine
    code_folder = code_folder_mac
    data_root = data_root_mac
end
data_root_children = 'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Child_MEG\'
data_root_adults = 'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\'




%% Add in necessary code to path
addpath(fullfile(code_folder,"\Matlab_toolboxes\fieldtrip-master\fieldtrip-master2"))
ft_defaults
addpath(fullfile(code_folder,"\Matlab_toolboxes\fieldtrip-master\fieldtrip-master\external\eeglab"))
addpath(fullfile(code_folder,"\Matlab_toolboxes\fieldtrip-master\fieldtrip-master\external\mffmatlabio\private"))
addpath(genpath('E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Paul, Jordan\Jordan\Jordan Code\\'))

addpath(fullfile(code_folder,"\SequentialBayesianLearning-master\Sequential_Bayesian_learner-master"))
addpath(fullfile(code_folder,"\spm12-master\spm12-master"))

addpath(fullfile(code_folder,"\Matlab_toolboxes\eeglab\functions\popfunc"))
addpath(fullfile(code_folder,"\Matlab_toolboxes\eeglab\functions\adminfunc"))
addpath(fullfile(code_folder,"\Matlab_toolboxes\eeglab\functions\sigprocfunc"))
addpath(fullfile(code_folder,"\Matlab_toolboxes\eeglab\plugins\firfilt2.4"))

addpath(fullfile(code_folder,"\Roving_MMN_2020-master\Roving_MMN_2020-master\preprocessing_scripts_ME125_phase2_YSun"))
addpath(fullfile(code_folder,"\Matlab_toolboxes\unfold"))
addpath(genpath(fullfile(code_folder,'\Matlab_toolboxes\MVPA-Light-master')))
%%





%% Get all subjects
%% The script will automatically get the subjects from the current folder.
cd(matlab_general_code)
subject_rethm = [""];

folder_list = workspace_find_variable_list_files(data_root_adults);
for i=1:length(folder_list)
    folder = char(folder_list(i));
    regex = regexp(folder, '(\d\d\d\d)', 'match');
    if length(regex) > 0
        int_ = char(regex(1));
        subject_rethm(i) = int_;
    end
end
subject_rethm = unique(subject_rethm);
subject_rethm = rmmissing(subject_rethm); 

if length(who('-regexp','SUBJECTS_TO_RUN')) > 0
    subjects = [SUBJECTS_TO_RUN]
end

HGF346 = [];
HGF692 = [];
HGF1038 = [];
HGF1384 = [];
HGF1730 = [];
HGF2000 = [];
CS2000 = [];
BS2000 = [];
PS2000 = [];
issues_adults = [];
issues_sound_delay_adults = []
%% Do HGF for the same stimulus set that was given to all participants with code < 9000 (i.e. not done by Lance)
for j=1:length(subject_rethm)
    try
        
        subject = subject_rethm(j)
        folder = char(strcat(data_root_adults,subject, '\'))
        cd(folder)

        issue = 0;
%         hgf_346       = dir('*HGF346*.txt');
%         hgf_692       = dir('*HGF692*.txt');
%         hgf_1038      = dir('*HGF1038*.txt');
%         hgf_1384      = dir('*HGF1384*.txt');
%         hgf_1730      = dir('*HGF1730*.txt');
%         hgf_2000      = dir('*HGF2000*.txt');
        ps_2000      = dir('*PS2000.txt');
        cs_2000      = dir('*CS2000.txt');
        bs_2000      = dir('*BS2000.txt');

        hgf_346       = dir('*HGF346 PE2_new.txt');
        hgf_346       = dir('*HGF346 PE2_mod.txt');
        hgf_346       = dir('*HGF346 PE2_mod_baked v2.txt');  

        hgf_692       = dir('*HGF692 PE2_new.txt');
        hgf_692       = dir('*HGF692 PE2_mod.txt');
        hgf_692       = dir('*HGF692 PE2_mod_baked v2.txt');

        hgf_1038      = dir('*HGF1038 PE2_new.txt');
        hgf_1038      = dir('*HGF1038 PE2_mod.txt');    
        hgf_1038      = dir('*HGF1038 PE2_mod_baked v2.txt');

        hgf_1384      = dir('*HGF1384 PE2_new.txt');
        hgf_1384      = dir('*HGF1384 PE2_mod.txt');
        hgf_1384      = dir('*HGF1384 PE2_mod_baked v2.txt');

        hgf_1730      = dir('*HGF1730 PE2_new.txt');   
        hgf_1730      = dir('*HGF1730 PE2_mod.txt');
        hgf_1730      = dir('*HGF1730 PE2_mod_baked v2.txt');

        hgf_2000      = dir('*HGF2000 PE2_new.txt');   
        hgf_2000      = dir('*HGF2000 PE2_mod.txt'); 
        hgf_2000      = dir('*HGF2000 PE2_mod_baked v2.txt');  

        sound_delay    = dir('*Sound Delay.txt');
        if length(sound_delay) < 1
            issues_sound_delay_adults = [issues_sound_delay_adults, subject]
        end
        if length(hgf_346) < 1
            HGF346 = [HGF346,subject];
            issue = issue + 1;
        end
        if length(hgf_692) < 1
            HGF692 = [HGF692,subject];
            issue = issue + 1;
        end
        if length(hgf_1038) < 1
            HGF1038 = [HGF1038,subject];
            issue = issue + 1
        end
        if length(hgf_1384) < 1
            HGF1384 = [HGF1384,subject];
            issue = issue + 1;
        end
        if length(hgf_1730) < 1
            HGF1730 = [HGF1730,subject];
            issue = issue + 1;
        end
        if length(hgf_2000) < 1
            HGF2000 = [HGF2000,subject];
            issue = issue + 1;
        end


%         if length(cs_2000) < 1
%             CS2000 = [CS2000,subject];
%             issue = issue + 1;
%         end
% 
%         if length(bs_2000) < 1
%             BS2000 = [BS2000,subject];
%             issue = issue + 1;
%         end
%         if length(ps_2000) < 1
%             PS2000 = [PS2000,subject];
%             issue = issue + 1;
%         end



        disp(subject)
        disp(hgf_346(1))
        disp(hgf_692(1))
        disp(hgf_1038(1))   
        disp(hgf_1384(1))
        disp(hgf_1730(1))
        disp(hgf_2000(1))

        
        tst = hgf_346(1).bytes >= hgf_692(1).bytes;
        tst = tst | hgf_692(1).bytes >= hgf_1038(1).bytes;
        tst = tst | hgf_1038(1).bytes >= hgf_1384(1).bytes;
        tst = tst | hgf_1384(1).bytes >= hgf_1730(1).bytes; % Note, in some cases there are false positive 'issues' when the sequence recorded was shorter and hence e.g. HGF1384 is the same size as HGF1730
        tst = tst | hgf_1730(1).bytes > hgf_2000(1).bytes;
        issue = issue + tst;
       



    catch
        %disp("Issue")
        issues_adults = [issues_adults, subject];
        issue = 1;
    end
    if issue > 0
        issues_adults = [issues_adults, subject];
    end
end
issues_adults = unique(issues_adults)




%% Get all children
%% The script will automatically get the subjects from the current folder.
cd(matlab_general_code)
subject_rethm = [""];

folder_list = workspace_find_variable_list_files(data_root_children);
for i=1:length(folder_list)
    folder = char(folder_list(i));
    regex = regexp(folder, '(\d\d\d\d)', 'match');
    if length(regex) > 0
        int_ = char(regex(1));
        subject_rethm(i) = int_;
    end
end
subject_rethm = unique(subject_rethm);
subject_rethm = rmmissing(subject_rethm); 

HGF346_children = [];
HGF692_children = [];
HGF1038_children = [];
HGF1384_children = [];
HGF1730_children = [];
HGF2000_children = [];
CS2000_children = [];
BS2000_children = [];
PS2000_children = [];
issues_children = [];
issues_sound_delay_kids = []
%% Do HGF for the same stimulus set that was given to all participants with code < 9000 (i.e. not done by Lance)
for j=1:length(subject_rethm)
    try
        subject = subject_rethm(j)
        folder = char(strcat(data_root_children,subject, '\'))
        cd(folder)

        issue = 0;
%         hgf_346       = dir('*HGF346*.txt');
%         hgf_692       = dir('*HGF692*.txt');
%         hgf_1038      = dir('*HGF1038*.txt');
%         hgf_1384      = dir('*HGF1384*.txt');
%         hgf_1730      = dir('*HGF1730*.txt');  
%         hgf_2000     = dir('*HGF2000*.txt');  
        ps_2000      = dir('*PS2000.txt');
        cs_2000      = dir('*CS2000.txt');
        bs_2000      = dir('*BS2000.txt');

        hgf_346       = dir('*HGF346 PE2_new.txt');
        hgf_346       = dir('*HGF346 PE2_mod.txt');
        hgf_346       = dir('*HGF346 PE2_mod_baked v2.txt');  

        hgf_692       = dir('*HGF692 PE2_new.txt');
        hgf_692       = dir('*HGF692 PE2_mod.txt');
        hgf_692       = dir('*HGF692 PE2_mod_baked v2.txt');

        hgf_1038      = dir('*HGF1038 PE2_new.txt');
        hgf_1038      = dir('*HGF1038 PE2_mod.txt');    
        hgf_1038      = dir('*HGF1038 PE2_mod_baked v2.txt');

        hgf_1384      = dir('*HGF1384 PE2_new.txt');
        hgf_1384      = dir('*HGF1384 PE2_mod.txt');
        hgf_1384      = dir('*HGF1384 PE2_mod_baked v2.txt');

        hgf_1730      = dir('*HGF1730 PE2_new.txt');   
        hgf_1730      = dir('*HGF1730 PE2_mod.txt');
        hgf_1730      = dir('*HGF1730 PE2_mod_baked v2.txt');

        hgf_2000      = dir('*HGF2000 PE2_new.txt');   
        hgf_2000      = dir('*HGF2000 PE2_mod.txt'); 
        hgf_2000      = dir('*HGF2000 PE2_mod_baked v2.txt');  
     
        sound_delay    = dir('*Sound Delay.txt');
        if length(sound_delay) < 1
            issues_sound_delay_kids = [issues_sound_delay_kids, subject]
        end
        if length(hgf_346) < 1
            HGF346_children = [HGF346_children,subject];
            issue = issue + 1;
        end
        if length(hgf_692) < 1
            HGF692_children  = [HGF692_children,subject];
            issue = issue + 1;
        end
        if length(hgf_1038) < 1
            HGF1038_children  = [HGF1038_children,subject];
            issue = issue + 1;
        end
        if length(hgf_1384) < 1
            HGF1384_children  = [HGF1384_children,subject];
            issue = issue + 1;
        end
        if length(hgf_1730) < 1
            HGF1730_children  = [HGF1730_children,subject];
            issue = issue + 1;
        end
        if length(hgf_2000) < 1
            HGF2000_children  = [HGF2000_children,subject];
            issue = issue + 1;
        end

        if length(cs_2000) < 1
            CS2000_children  = [CS2000_children,subject];
            issue = issue + 1;
        end

        if length(bs_2000) < 1
            BS2000_children  = [BS2000_children,subject];
            issue = issue + 1;
        end
        if length(ps_2000) < 1
            PS2000_children  = [PS2000_children,subject];
            issue = issue + 1;
        end

        tst = hgf_346(1).bytes >= hgf_692(1).bytes;
        tst = tst | hgf_692(1).bytes >= hgf_1038(1).bytes;
        tst = tst | hgf_1038(1).bytes >= hgf_1384(1).bytes;
        tst = tst | hgf_1384(1).bytes >= hgf_1730(1).bytes; % Note, in some cases there are false positive 'issues' when the sequence recorded was shorter and hence e.g. HGF1384 is the same size as HGF1730
        tst = tst | hgf_1730(1).bytes > hgf_2000(1).bytes;
        issue = issue + tst;



    catch
        %disp("Issue")
        issue = 1;
        issues_children = [issues_children, subject];
    end
    if issue > 0
        issues_children = [issues_children, subject];
    end
end
issues_children = unique(issues_children)
