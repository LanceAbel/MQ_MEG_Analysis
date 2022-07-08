SUBJECTS = {'2696','2697'}% '2678'}


% Generates HGF, BS, CS, PS predictors
% After running this, fix output with compare_timing_v3.py
NUM_TONES_LONG = 692 % Length of predictor array to generate



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
data_root = char(tlines(4));
data_root_mac = char(tlines(5));
matlab_general_code = char(tlines(6));
code_folder = char(tlines(7));
code_folder_mac = char(tlines(8));
if ~windows_machine
    code_folder = code_folder_mac
    data_root = data_root_mac
end




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
% Folders to look in:
%folder_list = workspace_find_variable_list_files('E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Child_MEG\\');
folder_list = workspace_find_variable_list_files(data_root);

for i=1:length(folder_list)
    folder = char(folder_list(i));
    regex = regexp(folder, '(\d\d\d\d)', 'match');
    if length(regex) > 0
        int_ = char(regex(1))
        subject_rethm(i) = int_;
    end
end
subject_rethm = unique(subject_rethm);
subject_rethm = rmmissing(subject_rethm) 
if length(who('-regexp','SUBJECTS')) > 0
    subject_rethm = SUBJECTS
end
%
% new_subjects = [""] % Sets it as string
% for j=1:length(subject_rethm)
%     new_subjects(j) = arrayfun(@num2str,subject_rethm(j),'un',0) 
% end
%new_subjects = rmmissing(new_subjects) 
%subject_rethm = new_subjects



issues_with_hgf = []
issues_with_binary = []
stim_folder = 'E:\oDrive\OneDrive\Docs Sync\Jobs And Money\Careers\Study\Uni\MAC\MQ MRes\Coding\Experiment\New Experiment\Tone Inputs\\'
cd(stim_folder)
all_stim_sequences = dir('*.txt');
txtfiles = all_stim_sequences    % Gets first one somehow

%% Do HGF for the identical stimulus set that was given to all participants with code 9000+


% for j = 3:3
% for j=1:length(txtfiles)
% 
%     clear PS1 CS1 BS1 PE2 PE3 PWPE2 PWPE3 stim stim_before
%     txtfile = txtfiles(j)
%     txtfile_name = fullfile(txtfile.folder,txtfile.name);
%     TF = contains(txtfile_name, 'HGF'); % Don't run on generated files
%     if TF==0
%         fid=fopen(txtfile_name);
%         tline = fgetl(fid);
%         freqs = cell(0,1);
%         while ischar(tline)
%             freqs{end+1,1} = tline;
%             tline = fgetl(fid);
%         end
%         freqs = cellfun(@str2num,freqs)
%         fclose(fid);
%     
%     
%         % The HGF ignores actual tone frequencies, you must input 1, 2... N where there are N tones
%         sorted_freqs = sort(freqs);
%         sorted_freqs = unique(sorted_freqs)';
%         tones_categoricals = 1:length(sorted_freqs); 
%     
%      
%         % Write tone identities
%         identities_filename = strcat(txtfile.name, ' identities')
%         identities_filename = strrep(identities_filename, '.txt', '')
%         identities_filename = strcat(identities_filename, ' .txt')
%         tones_train = [];
%         txtfile_name_new = fullfile(txtfile.folder, identities_filename);
%         %fileID = fopen(txtfile_name_new, 'w');    
%         for r=1:NUM_TONES_LONG % length(freqs)
%             freq_to_match = freqs(r);
%             index = find(sorted_freqs==freq_to_match);
%             tones_categorical = fix(tones_categoricals(index));
%             tones_categorical = int64(tones_categorical);
%             tones_train = [tones_train num2str(tones_categorical) newline];
%             %tones_train = strcat(tones_train,num2str(tones_categorical), '\n');
%             %fprintf(fileID,'%.0f\n', tones_categorical);
%             %fclose(fileID);
%         end
%     
%         % Write file corresponding to identities
%         %writematrix(tones_train, identities_filename);
%         tones_train_split = strsplit(tones_train, newline) 
%         tones_train_split(length(tones_train_split)) = [] % Remove last
%         to_fit_to_hgf = str2double(tones_train_split)
%         to_fit_to_hgf = to_fit_to_hgf'                  % Transpose
%         % Fit HGF
%         %est_tones_hgf_whatworld_mod_baked = tapas_fitModel([], to_fit_to_hgf, 'tapas_hgf_whatworld_config', 'tapas_bayes_optimal_whatworld_config');
%         est_tones_hgf_whatworld_mod_baked = tapas_fitModel([], to_fit_to_hgf, 'tapas_hgf_whatworld_config_mod', 'tapas_bayes_optimal_whatworld_config')        
%         % Write HGF output
%     
%         HGF_name = ' HGF_Settings999_ '
%         SeqLength = num2str(NUM_TONES_LONG)
%         Ext = '.mat'
%         hgf_savename = strcat(txtfile.name, HGF_name, SeqLength)
%         hgf_savename = strrep(hgf_savename, '.txt', '.mat')
%         hgf_savename = strrep(hgf_savename, '.mat', '')
%         hgf_savename = strcat(hgf_savename, '.mat')
%     
%         save(hgf_savename,  'est_tones_hgf_whatworld_mod_baked');  
%         for i = 2:length(est_tones_hgf_whatworld_mod_baked.u_orig)
%            stim_before  = est_tones_hgf_whatworld_mod_baked.u_orig(i-1);
%            stim         = est_tones_hgf_whatworld_mod_baked.u_orig(i);
%            PE2(i)       = (est_tones_hgf_whatworld_mod_baked.traj.da(i-1,1,stim,stim_before));    
%            PE3(i)       = sum(est_tones_hgf_whatworld_mod_baked.traj.da(i-1,2,stim_before,:),'omitnan');
%            PWPE2(i) = (est_tones_hgf_whatworld_mod_baked.traj.epsi(i-1,2,stim,stim_before));       
%            PWPE3(i) = sum(est_tones_hgf_whatworld_mod_baked.traj.epsi(i-1,3,stim_before,:),'omitnan');   
%         end
%     
%         newtxtfile_name = strrep(hgf_savename, Ext, '');
%         newtxtfile_name = strcat(newtxtfile_name, ' PE2.txt');
%         writematrix(PE2, newtxtfile_name);
%         newtxtfile_name = strrep(hgf_savename, Ext, '');
%         newtxtfile_name = strcat(newtxtfile_name, ' PE3.txt');
%         writematrix(PE3, newtxtfile_name);
%         newtxtfile_name = strrep(hgf_savename, Ext, '');
%         newtxtfile_name = strcat(newtxtfile_name, ' PWPE2.txt');
%         writematrix(PWPE2, newtxtfile_name);
%         newtxtfile_name = strrep(hgf_savename, Ext, '');
%         newtxtfile_name = strcat(newtxtfile_name, ' PWPE3.txt');
%         writematrix(PWPE3, newtxtfile_name);
%     end
% end



%% Do HGF for the distinct stimulus sets that were given to each participant with code < 9000 (i.e. not done by Lance)
for j=1:length(subject_rethm)
    try
        clear PS1 CS1 BS1 PE2 PE3 PWPE2 PWPE3 stim stim_before
        subject = subject_rethm(j)
        disp(subject)
        folder = char(strcat(data_root,subject, '\'))
        cd(folder)

        stimfile_norm       = dir('*oddball*run*.txt');
        stimfile_alt        = dir('*event_meg*.txt');
        % CAN NOT RELY ON THE ODDBALL_RUN FILE, SO SET to use stimfile_alt.
        stimfile_norm = []
        if length(stimfile_norm) == 0
            stimfile        = stimfile_alt;
        else
            stimfile        = stimfile_norm;
            tones           = load(fullfile(stimfile.folder,stimfile.name));
            tones_category  = tones(:,[3]); % Gets the third column, which contains the times
        end
        
        
        if length(stimfile_norm) == 0 % Then we need to examine only the first line
            fid=fopen(fullfile(stimfile.folder,stimfile.name));
            tline = fgetl(fid);
            tones = cell(0,1);
            while ischar(tline)
                try
                    tmp = strsplit(tline, ',');
                    tone_cell = tmp(1,3);
                    tone_val = cellfun(@str2num,tone_cell);
                    tones{end+1,1} = tone_val;
                catch
                    disp("Skipping first line");
                end
                tline = fgetl(fid);
            end
            fclose(fid);
            tones(1) = [];
            tones = tones(1:min(length(tones),NUM_TONES_LONG));
        else
            tones = tones_category(1:min(length(tones_category),NUM_TONES_LONG));
        end
        
        if length(stimfile_norm) == 0
            tones = cell2mat(tones);
        end
        
        % Fit HGF
        est_tones_hgf_whatworld_mod_baked = tapas_fitModel([], tones, 'tapas_hgf_whatworld_config_mod', 'tapas_bayes_optimal_whatworld_config')
        % Write HGF output
        HGF_name = 'HGF';
        SeqLength = num2str(NUM_TONES_LONG);
        Ext = '.mat';
        hgf_savename = strcat(HGF_name, SeqLength); % HGF346.txt
        hgf_savename = strrep(hgf_savename, '.txt', '.mat'); % HGF346.mat
        hgf_savename = strrep(hgf_savename, '.mat', ''); %HGF346
        hgf_savename = strcat(hgf_savename, '.mat'); % HGF346.mat
    
        hgf_savename2 = strrep(hgf_savename, '.mat', '_mod_baked v2.mat');
        %save(hgf_savename,  'est_tones_hgf_whatworld_mod_baked');  
        save(hgf_savename2,  'est_tones_hgf_whatworld_mod_baked');  

        for i = 2:length(est_tones_hgf_whatworld_mod_baked.u_orig)

           stim_before  = est_tones_hgf_whatworld_mod_baked.u_orig(i-1);
           stim         = est_tones_hgf_whatworld_mod_baked.u_orig(i);
           PE2(i)       = (est_tones_hgf_whatworld_mod_baked.traj.da(i-1,1,stim,stim_before));
           PWPE2(i)     = (est_tones_hgf_whatworld_mod_baked.traj.epsi(i-1,2,stim,stim_before));  
           PE3(i)       = sum(est_tones_hgf_whatworld_mod_baked.traj.da(i-1,2,stim_before,:),'omitnan');
           PWPE3(i)     = sum(est_tones_hgf_whatworld_mod_baked.traj.epsi(i-1,3,stim_before,:),'omitnan');   
        end
    
        newtxtfile_name = strrep(hgf_savename, Ext, '');
        newtxtfile_name = strcat(newtxtfile_name, ' PE2_mod_baked v2.txt');
        PE2 = PE2/max(PE2) % Single value has to be overwritten
        writematrix(PE2, newtxtfile_name);
        newtxtfile_name = strrep(hgf_savename, Ext, '');
        newtxtfile_name = strcat(newtxtfile_name, ' PE3_mod_baked v2.txt');
        writematrix(PE3, newtxtfile_name);
        newtxtfile_name = strrep(hgf_savename, Ext, '');
        newtxtfile_name = strcat(newtxtfile_name, ' PWPE2_mod_baked v2.txt');
        PWPE2 = PWPE2/max(PWPE2) % Single value has to be overwritten
        writematrix(PWPE2, newtxtfile_name);
        newtxtfile_name = strrep(hgf_savename, Ext, '');
        newtxtfile_name = strcat(newtxtfile_name, ' PWPE3_mod_baked v2.txt');    
        writematrix(PWPE3, newtxtfile_name);



%% Do for cs/ps/bs
%
%         if NUM_TONES_LONG == 2000
%             try
%                 %% Fit binary methods based on whether it was a repeat or not
% 
%                 stimfile_norm       = dir('*oddball*run*.txt');  
%                 stimfile_alt        = dir('*event_sounds.txt');
%                 % CAN NOT RELY ON THE ODDBALL FILE, SO SET to use stimfile_alt.
%                 stimfile_norm = []
%                 if length(stimfile_norm) == 0
%                     stimfile        = stimfile_alt;
%                 else
%                     stimfile        = stimfile_norm;
%                     tones           = load(fullfile(stimfile.folder,stimfile.name));
%                     tones_repeats  = tones(:,[4]); % Gets the third column, which contains the times
%                 end
%                 
%                 
%                 if length(stimfile_norm) == 0 % Then we need to examine excluding the first line
%                     fid=fopen(fullfile(stimfile.folder,stimfile.name));
%                     tline = fgetl(fid);
%                     tones = cell(0,1);
%                     while ischar(tline)
%                         try
%                             tmp = strsplit(tline, ',');
%                             tone_cell = tmp(1,3);
%                             tone_val = cellfun(@str2num,tone_cell);
%                             tones{end+1,1} = tone_val;
%                         catch
%                             disp("Skipping first line");
%                         end
%                         tline = fgetl(fid);
%                     end
%                     fclose(fid);
%                     tones(1) = [];
%                     tones = tones(1:min(length(tones),NUM_TONES_LONG));
%                 else
%                     tones = tones_repeats(1:min(length(tones_repeats),NUM_TONES_LONG));
%                 end
%                 
%                 if length(stimfile_norm) == 0
%                     tones = cell2mat(tones);
%                 end
%         
%         
%         
%                 seq  = [tones]==1; %find sequence starts
%                 seq2 = [];
%                 
%                 %Binary toggle between states
%                 d = 1;
%                 for i=1:length(seq)
%                     if seq(i)==1
%                         d = d+1;
%                     end
%                     seq2(i) = mod(d,2);
%                 end
%     
%                 
%                 % Calculate three different types of surprisal based on the binarised input: Predictive surprisal, Bayesian Surprisal, Confidence-corrected surprisal
%                 [PS1,BS1,CS1] = BL_Betabern_TP(seq2, 0.5);
%                 %https://github.com/kathrintertel/Sequential_Bayesian_learner
%     
%                 PredictorName = 'PS ';
%                 Ext = num2str(NUM_TONES_LONG)
%                 newtxtfile_name = strcat(PredictorName, Ext, '.txt');
%                 writematrix(PS1, newtxtfile_name);
%                 PredictorName = 'BS ';
%                 Ext = num2str(NUM_TONES_LONG)
%                 newtxtfile_name = strcat(PredictorName, Ext, '.txt');
%                 writematrix(BS1, newtxtfile_name);
%                 PredictorName = 'CS ';
%                 Ext = num2str(NUM_TONES_LONG)
%                 newtxtfile_name = strcat(PredictorName, Ext, '.txt');
%                 writematrix(CS1, newtxtfile_name);
%     
%     
%     
%     
%             catch
%                 issues_with_binary = [issues_with_binary, convertCharsToStrings(subject)] 
%             end
%         end

    catch
        issues_with_hgf = [issues_with_hgf, convertCharsToStrings(subject)]    
    end        

end    










