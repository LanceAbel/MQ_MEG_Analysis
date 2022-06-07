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
subject = char(tlines(1));
type = char(tlines(2));
folder_data = char(tlines(3));
data_root = char(tlines(4));
data_root_mac = char(tlines(5));
matlab_general_code = char(tlines(6));
code_folder = char(tlines(7));
code_folder_mac = char(tlines(8));
experiment_number = char(tlines(9));
num_channels_adult =  str2double(tlines(10));
num_channels_child =  str2double(tlines(11));
audio_channel_adult =  str2double(tlines(12));
audio_channel_child =  str2double(tlines(13));
first_trigger_adult =  str2double(tlines(14));
first_trigger_child =  str2double(tlines(15));
NUM_TONES =  str2double(tlines(16));
TRIAL_NUMBERS = 1:NUM_TONES;
if ismember(type,['Lance','Adult'])
    adult = 1;
    load(fullfile(data_root,'neighbours_160.mat'));
    num_channels = num_channels_adult;
else
    adult = 0;
    load(strcat(data_root,'neighbours_125.mat'));
    num_channels = num_channels_child;
end

if ~windows_machine
    code_folder = code_folder_mac
end
folder = strcat(data_root,subject) % folder with subject subfolders
cd(matlab_general_code)
add_all_paths
cd(folder)

    

output_path = strcat(folder,'CHILD_MEMES_V2\')
folder_meg_research = strcat(code_folder, 'Macquarie-MEG-Research')
% Path to MRI Library for MEMES
path_to_MRI_library = strcat(folder_data,'\Child_MRI\new_HCP_library_for_MEMES\\') ;
% Path to MQ_MEG_Scripts
path_to_MQ_MEG_Scripts = strcat(folder_meg_research,'\MQ_MEG_Scripts-master\MQ_MEG_Scripts-master\') 
% Path to MEMES
path_to_MEMES = strcat(folder_meg_research,'\MEMES-master\MEMES-master\') 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mq_preprocessing_example.m 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Set up paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path to where the data should be saved (make sure this ends with /)
if windows_machine
    save_path   = strcat(folder);%, '\sample_data_processed\');  % '/Volumes/Robert T5/sample_data_processed/';
else
    save_path   = strcat(folder);%, '/sample_data_processed/');  % '/Volumes/Robert T5/sample_data_processed/';   
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Add MQ_MEG_Scripts and MEMES to path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Adding MQ_MEG_Scripts and MEMES to your MATLAB path');
warning(['Please note that MQ_MEG_Scripts and MEMES are designed for'...
    ' MATLAB 2016b or later and have been tested using Fieldtrip'...
    ' version 20181213']);
addpath(genpath(path_to_MQ_MEG_Scripts));
addpath(genpath(path_to_MEMES));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Specifiy Subject ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%subject = subject %  '3566'; % LANCE: DEFINED ABOVE

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Make a subject specific results folder for saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Making subject specific folder for saving');

% Get the path to the saving directory
dir_name = [save_path];
% Make the directory!
mkdir(dir_name);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Specify paths to confile, mrkfile, elpfile and hspfile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % LANCE: Need to automate this to work on any subject
%for i = 2:length(act_subject)
cd(folder)
pathname = folder

elpfile = dir('*.elp').name
hspfile = dir('*.hsp').name

all_mrks = dir('*.mrk')
all_mrks_table = struct2table(all_mrks); % convert the struct array to a table
all_mrks_by_size = sortrows(all_mrks_table, 'date', 'ascend'); % sort the table by 'size', as 15 min roving oddball size will be larger than 10 min resting state
all_mrks_by_size = table2struct(all_mrks_by_size);
mrkfile = all_mrks_by_size.name

% confile =   '2872_CM_ME125_2018_01_12_B1_denoise_rethm.con'
if adult == 1
    all_cons = dir(strcat('*experiment',experiment_number,'_tspca.con'));
else
    all_cons = dir(strcat('*experiment',experiment_number,'_rethm_tspca.con')); % dir('*experiment1_rethm_tspca.con');
    if length(all_cons) == 0
        all_cons = dir('*denoise_rethm.con');
    end
end

% if length(all_cons) == 0
%     if adult == 1
%         all_cons = dir('*B3.con');
%     elseif adult == 0
%         all_cons = dir('*B1.con');
%     end
% end
if length(all_cons) == 0
    all_cons = dir('*.con');
end
if length(all_cons) == 0
    disp("NO CON FILE FOUND")
end
all_cons_table = struct2table(all_cons); % convert the struct array to a table
all_cons_by_size = sortrows(all_cons_table, 'bytes', 'descend'); % sort the table by 'size', as 15 min roving oddball size will be larger than 10 min resting state
all_cons_by_size = table2struct(all_cons_by_size);
confile = all_cons_by_size.name    % Gets first one somehow



if adult == 0
    child_MEMES_v2(dir_name,elpfile,hspfile,confile,mrkfile,path_to_MRI_library)
end


% Get the path to the saving directory
dir_name = [save_path];
cd(dir_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Check for saturations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if adult
    [sat] = mq_detect_saturations(dir_name,confile,0.01,'adult','array'); % Sourced from HTML, is this more up to date?
else
    [sat] = mq_detect_saturations(dir_name,confile,0.01,'child','array'); % Sourced from HTML, is this more up to date?
end
save sat sat



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. Downsample Headshape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% LANCE - downsample_headshape NO LONGER WORKS,
%headshape_downsampled = downsample_headshape(hspfile,'yes',0); 
%headshape_downsampled = downsample_headshape(hspfile,include_facial_points='yes',remove_zlim=0,remove_below_nasion=0); 
% LANCE - try downsample_headshape_new
cfg = []
cfg.facial_info = 'yes';
cfg.remove_zlim = 0;
cfg.downsample_facial_info  = 'yes';
%cfg.method = 'nonuniform' % Lance: 'gridaverage' produces error Incorrect Method Specified Unrecognized function or variable 'decimated_headshape'.


%If your data has extensive facial information, please consider using the additional options to match with the MEMES database:
cfg.facial_info_above_z           = 20 % remove points Xmm above the nasion on the z-axis (up-down)
cfg.facial_info_below_z           = 20 % 80 % remove points Xmm below the nasion on the z-axis (up-down)
cfg.facial_info_above_y           = 70 % remove points Xmm above the nasion on the y-axis (left-right) 
cfg.facial_info_above_y           = 70 % remove points Xmm below the nasion on the y-axis (left-right) 
cfg.facial_info_below_x           = 20 % remove points Xmm below the nasion on the x-axis (forwards-backwards)
headshape_downsampled = downsample_headshape_new(cfg, hspfile); % LANCE - REPLACED BY
figure; ft_plot_headshape(headshape_downsampled);


% Save
cd(dir_name);
disp('Saving headshape_downsampled');
save headshape_downsampled headshape_downsampled;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Realign Sensors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[grad_trans] = mq_realign_sens(dir_name,elpfile,hspfile,...
    confile,mrkfile,'','rot3dfit');

print('grad_trans','-dpng','-r200');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9. Read & Plot ReTHM Data - for children only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if adult == 0
    cd(strcat(code_folder, '\Macquarie-MEG-Research\MQ_MEG_Scripts-master\MQ_MEG_Scripts-master\Tools\reTHM'))
    bad_coil = '';
    [head_movt, confound]  = get_reTHM_data(dir_name,confile,grad_trans,...
        hspfile,bad_coil,0.99)
    save head_movt head_movt
    save confound confound
end

disp("Enter in bad coils")
%pause
[grad_trans] = mq_realign_sens(dir_name,elpfile,hspfile,...
    confile,mrkfile,{'LPAred'},'icp');
%%




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Read in raw MEG data & apply filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('Running Preprocessing Script for Project ME176 - Alien Task');
%disp('Running Preprocessing Script for Project ME125 - Roving oddball task');
disp('Running Preprocessing Script for Project ME199 - Roving oddball task');

% CD to correct directory
cd(dir_name);

% Epoch the whole dataset into one continous dataset and apply
% the appropriate filters
cfg = [];
cfg.headerfile = confile;
cfg.datafile = confile;
cfg.trialdef.triallength = Inf;
cfg.trialdef.ntrials = 1;
cfg = ft_definetrial(cfg);

cfg.continuous = 'yes';
alldata = ft_preprocessing(cfg);

% Band-pass filter between 0.5-250Hz to help visualisation
cfg.continuous = 'yes';
cfg.bpfilter = 'yes';
cfg.bpfreq = [0.5 250];
disp("Bandpass Filter")
alldata = ft_preprocessing(cfg);

% Deal with 50Hz line noise using a bandstop filter
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = [49.5 50.5];
disp("Dealing with line noise 50Hz")
alldata = ft_preprocessing(cfg,alldata);

% Deal with 100Hz line noise using a bandstop filter
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = [99.5 100.5];
disp("Dealing with line noise 100Hz")
alldata = ft_preprocessing(cfg,alldata);

ft_databrowser(cfg, alldata);


% Create layout file for later and save
cfg             = [];
cfg.grad        = alldata.grad;
lay             = ft_prepare_layout(cfg, alldata);
save lay lay

% Cut out MEG channels
cfg = [];
cfg.channel     = alldata.label(1:num_channels);
alldata         = ft_selectdata(cfg,alldata);


%% Lance cut out saturations: this caused errors on other script
if size(sat)>0 
satLabel=sat.label;
satTime=sat.time;
        for i=1:length(satLabel)
        [lia,locb]  = ismember(satTime{i},alldata.time{1,1});
        alldata.trial{1,1}(:,locb)=NaN;
        end 
end
%%




L = length(alldata.trial{1,1}(:,1)); % length(data.label

%% LANCE - Matlab artefact removal settings
%% LANCE - show times of likely muscle artefacts. Run before downsampling given filter freq
%The following settings are useful for identifying muscle artifacts:
muscledata = alldata;
cfg = [];
cfg.preproc.bpfilter    = 'yes';
cfg.preproc.bpfreq      = [110 140];
cfg.preproc.bpfiltord   =  8;
cfg.preproc.bpfilttype  = 'but';
cfg.preproc.rectify     = 'yes';
cfg.preproc.boxcar      = 0.2;
muscle = ft_preprocessing(cfg,muscledata);
ft_databrowser(cfg, muscle);


% Lance - custom artefact removal
outlier_z = 4; % what z score to look for as possible outlier
min_consecutive = 10; % How many data points in a row must exceed (positive or negative) z-score in order to label it an artefact
min_time_between_artefacts = 5; % How many seconds since the last artefact must have passed before we label it a "new" artefact
data_blocks = cell2mat(muscle.trial(1));
outliers_muscle = struct('overshoots_muscle', num2cell(1:L), ...
            'undershoots_muscle', num2cell(1:L), ...
            'undershoots_muscle_indices', num2cell(1:L), ...
            'overshoots_muscle_indices',num2cell(1:L));

indices_undershoot_muscle = struct('undershoots_muscle', num2cell(1:L));
indices_overshoot_muscle = struct('overshoots_muscle', num2cell(1:L));
for chan = 1:L
    %disp(chan)    
    data_chan   = data_blocks(chan,:);
    %disp(min(data_chan))
    %disp(max(data_chan))    
    avg         = nanmean(data_chan);
    stdev       = nanstd(data_chan);
    lower_bound = avg-outlier_z*stdev;
    upper_bound = avg+outlier_z*stdev;
   
    small_values = data_chan(data_chan < lower_bound);
    large_values = data_chan([data_chan] > upper_bound);

    %Bools
    %For this channel
    undershoots_muscle = data_chan(1,:) < lower_bound;
    overshoots_muscle = data_chan(1,:) > upper_bound;
    %sum(undershoots_muscle)
    %sum(overshoots_muscle)
    %Build up list for every channel
    outliers_muscle(chan).undershoots_muscle = undershoots_muscle;
    outliers_muscle(chan).overshoots_muscle = overshoots_muscle;
    outliers_muscle(chan).undershoots_muscle_indices = find(undershoots_muscle)
    outliers_muscle(chan).overshoots_muscle_indices = find(overshoots_muscle)
    outliers_muscle(chan).undershoots_muscle_times = [];
    outliers_muscle(chan).overshoots_muscle_times = []; 
    for position = 1:length(outliers_muscle(chan).undershoots_muscle_indices)
        index = outliers_muscle(chan).undershoots_muscle_indices(position);
        time = alldata.time{1,1}(index);
        outliers_muscle(chan).undershoots_muscle_times = [outliers_muscle(chan).undershoots_muscle_times, time];
        outliers_muscle(chan).overshoots_muscle_times = [outliers_muscle(chan).overshoots_muscle_times, time];       
    end
    % Get start time of overshoot/undershoot
    outliers_muscle(chan).undershoots_muscle_unique_times = []
    if length(outliers_muscle(chan).undershoots_muscle_times) > 1
        outliers_muscle(chan).undershoots_muscle_unique_times = outliers_muscle(chan).undershoots_muscle_times(1)
    end
    for time_index = 2:length(outliers_muscle(chan).undershoots_muscle_times)
        time_gap = outliers_muscle(chan).undershoots_muscle_times(time_index) - outliers_muscle(chan).undershoots_muscle_times(time_index-1);
        if time_gap > min_time_between_artefacts
            outliers_muscle(chan).undershoots_muscle_unique_times = [outliers_muscle(chan).undershoots_muscle_unique_times,outliers_muscle(chan).undershoots_muscle_times(time_index)];
        end
    end
end

%%
% calculate channel-to-channel correlations
% for chan_outer =1:length(alldata.trial{1,1}(:,1)))
%     for chan_inner = 1:length(alldata.trial{1,1}(:,1))
corrs = corrcoef(alldata.trial{1,1}');
cfg             = [];
cfg.grad        = alldata.grad;
lay             = ft_prepare_layout(cfg, alldata);

num_channels = min(num_channels,length(lay.pos));
distances = zeros([num_channels,num_channels]);
for chan_outer =1:num_channels
    for chan_inner =1:num_channels
        x_distance = lay.pos(chan_inner,1)-lay.pos(chan_outer,1);
        y_distance = lay.pos(chan_inner,2)-lay.pos(chan_outer,2);
        distance = sqrt(x_distance^2 + y_distance^2);
        distances(chan_outer,chan_inner) = distance;
    end
end

distances_array = [];
correlations_array = [];
for chan_outer=1:num_channels
    for chan_inner = chan_outer:num_channels
        distance = distances(chan_outer,chan_inner);
        distances_array = [distances_array,distance];
        corr = corrs(chan_outer,chan_inner);
        correlations_array = [correlations_array, corr];
    end
end
disp(corrcoef(distances_array,correlations_array))
scatter(distances_array,correlations_array)

    
%%%%%BEGIN LANCE ADD IN STEP #10 - ICA FROM MEG_CHILD_EXAMPLE.M
% THIS IS AS WE WERE GETTING ERRORS in data_epoched = ft_redefinetrial(cfg,data_clean)
% ...which required % data_clean
% ...which is output after ICA in MEG_CHILD_EXAMPLE
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. Perform ICA on the continuous data & Remove ECG/EOG artefacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downsample Data
alldata_orig    = alldata; %save the original CLEAN data for later use 
cfg             = []; 
cfg.channel     = alldata.label(1:num_channels);
cfg.resamplefs  = 150; %downsample frequency 
cfg.detrend     = 'no'; 
disp('Downsampling data');
alldata         = ft_resampledata(cfg, alldata);
alldata_copy    = alldata;

%% LANCE ARTEFACT REJECTION - USE FT TOOLBOX - CAN REMOVE
cfg = [];
cfg.method                          = 'distance'; % 'triangulation' % 'distance'  % The  distance  method simply draws a circle of certain size around each sensor-position.  The radius of the circle is defined by cfg.neighbourdist. 'template'
cfg.layout                          = 'lay';
cfg.trl                             = alldata.trial;
cfg.headerfile                      = confile;
cfg.datafile                        = confile;
neighbours                          = ft_prepare_neighbours(cfg);
cfg.neighbours                      = neighbours;
cfg.feedback                        = 'yes';
cfg.continuous                      = 'yes' ;
 
% % TEST - remove jumps
cfg.artfctdef.zvalue.channel        = alldata_orig.cfg.channel;
cfg.artfctdef.zvalue.channel        = 'MEG';
cfg.artfctdef.zvalue.cutoff         = 3;
cfg.artfctdef.zvalue.artpadding     = 0; % 0.02 % seconds Get rid of areas where the artefact is beginning to appear but is sub detection threshold
cfg.artfctdef.zvalue.trlpadding     = 0; % 0.5 % seconds Get rid of the whole trial when there's an artefact found inside it?
cfg.artfctdef.zvalue.fltpadding     = 0; % 0.1 % seconds Edge-of-filter effects

cfg.artfctdef.zvalue.medianfiltord  = 9;
cfg.artfctdef.zvalue.cumulative     = 'yes';
cfg.artfctdef.zvalue.medianfilter   = 'yes';
cfg.artfctdef.zvalue.absdiff        = 'yes';
cfg.artfctdef.zvalue.interactive    = 'yes';
%[cfg, artifact_jump]                = ft_artifact_zvalue(cfg);

% 
% TEST - remove EOG
cfg.artfctdef.zvalue.channel        = 'EOG';
cfg.artfctdef.eog.bpfilter          = 'yes';
cfg.artfctdef.eog.bpfilttype        = 'but';
cfg.artfctdef.eog.bpfreq            = [2 15];
cfg.artfctdef.eog.bpfiltord         = 4;
cfg.artfctdef.eog.hilbert           = 'yes';
cfg.artfctdef.eog.cutoff            = 4;
cfg.artfctdef.eog.trlpadding        = 0.5;
cfg.artfctdef.eog.fltpadding        = 0.1;
cfg.artfctdef.eog.artpadding        = 0.1;
%[cfg, artifact_EOG]                 = ft_artifact_eog(cfg);

% TEST - remove muscle movement
cfg.artfctdef.zvalue.channel        = 'MRT*';
cfg.artfctdef.muscle.bpfilter       = 'yes'
cfg.artfctdef.muscle.bpfilttype     = 'but'
cfg.artfctdef.muscle.bpfreq         = [110 140]
cfg.artfctdef.muscle.bpfiltord      = 8;
cfg.artfctdef.muscle.hilbert        = 'yes';
cfg.artfctdef.muscle.cutoff         = 4;
cfg.artfctdef.muscle.boxcar         = 0.2;
cfg.artfctdef.muscle.trlpadding     = 0.1;
cfg.artfctdef.muscle.fltpadding     = 0.1;
cfg.artfctdef.muscle.artpadding     = 0.1;
%[cfg, artifact_muscle]              = ft_artifact_muscle(cfg);


% Then mplement the particular artefact rejection you want to
%cfg.artfctdef.jump.artifact         = artifact_jump;
%cfg.artfctdef.eog.artifact          = artifact_EOG; %
%cfg.artfctdef.muscle.artifact       = artifact_muscle;
cfg.artfctdef.reject                = 'complete';  % Remove whole trial
alldata                             = ft_rejectartifact(cfg, alldata);
%Inspect new data
% ft_databrowser(cfg, alldata);

all_data_blocks                     = alldata.trial{1,1};
L = length(all_data_blocks(:,1)) % length(data.label


% Examine differences in each channel at each timepoint
diffs_per_time = struct();
diffs_per_time.times = [repelem(0,length(all_data_blocks))];
diffs_found = struct();
diffs_found.chans = [repelem(diffs_per_time(),L)];
tic
for time_ = 1:length(all_data_blocks) % 9 minutes for 922k time samples, can downsample
    %disp(times_)
    channel_col_all = alldata_copy.trial{1,1}(:,time_);
    channel_col_all_z_reject = alldata.trial{1,1}(:,time_);    
    for chan = 1:L
        val_ = channel_col_all(chan);
        val_z_reject = channel_col_all_z_reject(chan);
        diff_bw = val_z_reject - val_;
        diffs_found.chans(chan).times(time_) = diff_bw;
    end
end
toc

% Test whether there is any difference between the signal before and after artefact rejection
any_diff = 0
for r = 1:L
    any_diff = any_diff + sum(diffs_found.chans(r).times);
end
disp(any_diff)
% Plot differences
for chan = 1:L
    max_ = max(abs(diffs_found.chans(chan).times(1:length(all_data_blocks))));   
    if max_ > 0.00000000001
        disp(chan)
        disp(max_)
        figure;
        plot(diffs_found.chans(chan).times(1:length(all_data_blocks)));
    end
end




%FT_CHANNELREPAIR repairs bad or missing channels in the data by replacing them with the plain average of of all neighbours, by a weighted average of all neighbours, by an  interpolation based on a surface Laplacian, or by spherical spline interpolating (see Perrin et al., 1989).
%Use as
%[interp] = ft_channelrepair(cfg, data)
% where the input data corresponds to the output from FT_PREPROCESSING.
% arft.badchannel                     = input('write badchannels: ');
% cfg.badchannel                      = arft.badchannel
% data_fixed                          = ft_channelrepair(cfg,data);








%% LANCE Identify times when MANY channels look off
%% LANCE - show times of likely EOG artefacts
eogdata = alldata
%The following settings are useful for identifying EOG artifacts:
cfg = [];
cfg.preproc.bpfilter    = 'yes'
cfg.preproc.bpfilttype  = 'but'
cfg.preproc.bpfreq      = [1 15]
cfg.preproc.bpfiltord   = 4
cfg.preproc.rectify     = 'yes'

eog = ft_preprocessing(cfg,eogdata);
ft_databrowser(cfg, eog);

data_blocks = cell2mat(eog.trial(1));
L = length(data_blocks(:,1)); % length(data.label
outliers_eog = struct('overshoots_eog', num2cell(1:L), ...
            'undershoots_eog', num2cell(1:L));
indices_undershoot_eog = struct('undershoots_eog', num2cell(1:L));
indices_overshoot_eog = struct('overshoots_eog', num2cell(1:L));
%outliers = struct('undershoots', containers.Map, 'overshoots', containers.Map)
for chan = 1:L
    %disp(chan)    
    data_chan   = data_blocks(chan,:);
    %disp(min(data_chan))
    %disp(max(data_chan))
    avg         = nanmean(data_chan);
    stdev       = nanstd(data_chan);
    lower_bound = avg-outlier_z*stdev;
    upper_bound = avg+outlier_z*stdev;
   
    small_values = data_chan(data_chan < lower_bound);
    large_values = data_chan([data_chan] > upper_bound);

    %Bools
    %For this channel
    undershoots_eog = data_chan(1,:) < lower_bound;
    overshoots_eog = data_chan(1,:) > upper_bound;
    %sum(undershoots_eog)
    %sum(overshoots_eog)
    %Build up list for every channel
    outliers_eog(chan).undershoots_eog = undershoots_eog;
    outliers_eog(chan).overshoots_eog = overshoots_eog;
    outliers_eog(chan).undershoots_eog_indices = find(undershoots_eog);
    outliers_eog(chan).overshoots_eog_indices = find(overshoots_eog);
    outliers_eog(chan).undershoots_eog_times = [];
    outliers_eog(chan).overshoots_eog_times = [];
    for position = 1:length(outliers_eog(chan).undershoots_eog_indices) % Populate the times that the outliers are found (secs)
        index = outliers_eog(chan).undershoots_eog_indices(position);
        time = alldata.time{1,1}(index);
        outliers_eog(chan).undershoots_eog_times = [outliers_eog(chan).undershoots_eog_times, time];
        outliers_eog(chan).overshoots_eog_times = [outliers_eog(chan).overshoots_eog_times, time];    
    end

    % Get start time of overshoot/undershoot
    outliers_eog(chan).undershoots_eog_unique_times = [];
    if length(outliers_eog(chan).undershoots_eog_times) > 1
        outliers_eog(chan).undershoots_eog_unique_times = outliers_eog(chan).undershoots_eog_times(1);
    end
    for time_index = 2:length(outliers_eog(chan).undershoots_eog_times)
        time_gap = outliers_eog(chan).undershoots_eog_times(time_index) - outliers_eog(chan).undershoots_eog_times(time_index-1);
        if time_gap > 3
            outliers_eog(chan).undershoots_eog_unique_times = [outliers_eog(chan).undershoots_eog_unique_times,outliers_eog(chan).undershoots_eog_times(time_index)];
        end
    end
end


% Calc stats for each channel
chan_stats = struct('average', num2cell(1:L), ...
                        'stdev', num2cell(1:L));
L = length(alldata.label);
for chan = 1:L
    channel_row = all_data_blocks(chan,:);
    avg_ = nanmean(channel_row);
    stdev_ = nanstd(channel_row);
    chan_stats(chan).average = avg_;
    chan_stats(chan).stdev = stdev_;
end
% Identify how many channels have deviations from average
multiple_deviations = zeros( [length(all_data_blocks),1]);
min_channels_with_devs = 4;
tic
for times_=1:length(all_data_blocks) % 9 minutes for 922k time samples, can downsample
    num_channels_with_devs = 0;
    channel_col = all_data_blocks(:,times_);
    for chan = 1:L
        val = channel_col(chan);
        diff = val - chan_stats(chan).average;
        diff_z = diff/chan_stats(chan).stdev;
        if diff_z > outlier_z || diff_z < -outlier_z
            num_channels_with_devs = num_channels_with_devs + 1;
        end
    end
    if num_channels_with_devs > min_channels_with_devs
        multiple_deviations(times_) = 1;
    end
end
sum(multiple_deviations) % How many distinct time points show deviations in many channels
toc

% % Look for a bunch of consecutive outliers
% times_undershoot_eog = []
% times_overshoot_eog = []
% time_last_undershoot = outliers_eog(chan).undershoots_eog_times(1)
% time_last_overshoot = outliers_eog(chan).undershoots_eog_times(1)
% for i = 1:length(outliers_eog(chan).undershoots_eog_times);
%     if outliers_eog(chan).undershoots_eog_times(i) -
%         consecutive_undershoot = consecutive_undershoot + 1;
%         if consecutive_undershoot > min_consecutive;
%             indices_undershoot_eog(chan) = [indices_undershoot_eog(chan),i];
%             indices_overshoot_eog(chan) = [indices_overshoot_eog(chan),i];  
%         else
%         end
%     else
%         consecutive_undershoot = 0;
%     end;
% end;



















% 
% 
% 
% % Run ICA
% disp('About to run ICA using the Runica method')
% cfg                 = [];
% %cfg.method          = 'runica'; % doesn't work
% cfg.method          = 'fastica';
% cfg.feedback        = 'textbar';
% cfg.numcomponent    = 50;
% cfg.maxNumIterations = 200
% 
% %cfg.channel = 'all'; 

% 
% comp                = ft_componentanalysis(cfg, alldata);
% cfg.method          = 'mtmfft';
% cfg.taper      = 'dpss'
% freq = ft_freqanalysis(cfg, comp);
% % 


alldata_orig    = alldata; %save the original CLEAN data for later use 
% Run ICA
disp('About to run ICA using the fastica method')
cfg                 = [];
cfg.method          = 'fastica';
cfg.numcomponent    = 50;
cfg.maxNumIterations = 200
cfg.feedback        = 'textbar';
comp                = ft_componentanalysis(cfg, alldata);

% Display Components - change layout as needed
cfg             = []; 
cfg.compscale   = 'local';
cfg.viewmode    = 'component'; 
cfg.layout      = lay;
cfg.position    = [1 1 800 700];
cfg.ylim        = [-1.149e-11 1.149e-11]
%cfg.ylim        = [ -2.8015e-11  5.7606e-11 ];
cfg.blocksize   = 5;
ft_databrowser(cfg, comp);

ft_hastoolbox('brewermap',1);
colormap123     = colormap(flipud(brewermap(64,'RdBu')));

% Decompose the original data as it was prior to downsampling 
disp('Decomposing the original data as it was prior to downsampling...');
cfg           = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_orig     = ft_componentanalysis(cfg, alldata_orig);






% Display Components - change layout as needed
cfg             = []; 
cfg.compscale   = 'local';
cfg.viewmode    = 'component'; 
cfg.layout      = lay;
%cfg.position    = [1 1 800 700];
cfg.ylim        = [-1.149e-11 1.149e-11]
%cfg.ylim        = [ -2.8015e-11  5.7606e-11 ];  % Components 9 & 10, and many 11-20 were huge Couldn't see others until changing axis
cfg.blocksize   = 5;
ft_databrowser(cfg, comp);




ft_hastoolbox('brewermap',1);
colormap123     = colormap(flipud(brewermap(64,'RdBu')));

% Decompose the original data as it was prior to downsampling 
disp('Decomposing the original data as it was prior to downsampling...');
cfg           = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_orig     = ft_componentanalysis(cfg, alldata_orig);







%% The original data can now be reconstructed, excluding specified components
% This asks the user to specify the components to be removed
disp('Enter components in the form [1 2 3]')
comp2remove     = []% input('Which components would you like to remove?\n');
cfg             = [];
cfg.component   = [comp2remove]; %these are the components to be removed

% %LANCE ADD
% cfg.layout = 'CTF151.lay'; % specify the layout file that should be used for plotting
% cfg.viewmode = 'component';
%ft_databrowser(cfg, comp);
% %END LANCE ADD

data_clean      = ft_rejectcomponent(cfg, comp_orig,alldata_orig);

%ft_databrowser(cfg, data_clean); % LANCE ADD

%%
%%%%% END LANCE ADD IN STEP #10







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. Epoch the data into trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    %cfg.continuous              = 'yes';       % NOT IN HTML
    cfg                         = [];
    cfg.dataset                 = confile;
    cfg.trialdef.prestim        = 0.1;          % LANCE FIX pre-stimulus interval
    cfg.trialdef.poststim       = 0.4;          % LANCE FIX post-stimulus interva
    cfg.trialdef.eventtype      = 'trigger';
    cfg.trialdef.eventvalue     = TRIAL_NUMBERS; % Trigger numbers
    cfg.Fs                      = 1000;
    
    if adult == 1
        % Adult
        cfg.First_Channel           = first_trigger_adult; % First trigger channel  146
        cfg.Last_Channel            = first_trigger_adult+NUM_TONES; % Last trigger channel   158
        cfg.Audio_Channel           = audio_channel_adult; % 0 for none  135
    else
        % Child
        cfg.First_Channel           = first_trigger_child; % First trigger channel  146
        cfg.Last_Channel            = first_trigger_child+NUM_TONES; % Last trigger channel   158
        cfg.Audio_Channel           = audio_channel_child; % 0 for none  135
    end
    
    %cfg.fixed_offset            = 244;   % Default = []; 42ms for those without audio channel.
    cfg.trialfun                = 'FindTriggers_AudioCh'; %
    cfg                         = ft_definetrial(cfg);
    data_epoched                = ft_redefinetrial(cfg,data_clean);
catch
    %cfg.continuous              = 'yes';       % NOT IN HTML
    cfg                         = [];
    cfg.dataset                 = confile;
    cfg.trialdef.prestim        = 0.1;          % LANCE FIX pre-stimulus interval
    cfg.trialdef.poststim       = 0.4;          % LANCE FIX post-stimulus interva
    cfg.trialdef.eventtype      = 'trigger';
    cfg.trialdef.eventvalue     = TRIAL_NUMBERS; % Trigger numbers
    cfg.Fs                      = 1000;
    
    if adult == 1
        % Adult
        cfg.First_Channel           = first_trigger_adult; % First trigger channel  146
        cfg.Last_Channel            = first_trigger_adult+NUM_TONES; % Last trigger channel   158
        cfg.Audio_Channel           = audio_channel_adult; % 0 for none  135
    else
        % Child
        cfg.First_Channel           = first_trigger_child; % First trigger channel  146
        cfg.Last_Channel            = first_trigger_child+NUM_TONES; % Last trigger channel   158
        cfg.Audio_Channel           = audio_channel_child; % 0 for none  135
    end
    
    cfg.fixed_offset            = 244;   % Default = []; For those without audio channel. 244 ms post RME upgrade (audio device 14, 42ms before RME upgrade)
    cfg.trialfun                = 'FindTriggers_AudioCh'; %
    cfg                         = ft_definetrial(cfg);
    data_epoched                = ft_redefinetrial(cfg,data_clean);
end



% Now we want to select deviant and predeviant trials of interest
% Get sequence of tones
event = data_epoched.cfg.event;
ggg = [];
for i = 1:length(event)
    ggg(i,1) = event(i).value;
end

% Find the deviant trials
deviant_trials = find(ggg == 1);

% Remove first trial
deviant_trials(1) = [];

% Now found deviant trials in which previous sequence length is greater
% than 3
fff = ggg(deviant_trials-1);
deviant_trials2 = deviant_trials(find(fff > 3));

% Get trial before (predeviant)
predeviant_trials = deviant_trials2-1;
fprintf('Found %d deviant trials\n', length(deviant_trials2));
fprintf('Found %d predeviant trials\n', length(predeviant_trials));

% Select data
cfg = [];
cfg.trials = deviant_trials2;
deviant = ft_selectdata(cfg,data_epoched);
cfg.trials = predeviant_trials;
predeviant = ft_selectdata(cfg,data_epoched);









%% START ADDITIONS FOR EXTRA DETAIL

%Understand triggers more
plot([event.sample], [event.value], '.')
%In case the triggers are of different types, e.g. like here
% disp(unique({event.type})) % 'trigger' only

%%
% search for "trigger" events
value  = [event(find(strcmp('trigger', {event.type}))).value]';
sample = [event(find(strcmp('trigger', {event.type}))).sample]';
% disp(value)  % Roving oddball so 1,2,3,4,5,1,2,3,4,1, ... etc
% disp(sample) % The times, in hundreds of seconds (start ~20sec, then every ~500msec)


% determine the number of samples before and after the trigger
cfg.trialdef.prestim        = 0.1;          % LANCE FIX pre-stimulus interval
cfg.trialdef.poststim       = 0.4;          % LANCE FIX post-stimulus interva
pretrig  = -round(cfg.trialdef.prestim * 1000); % pretrig  = -round(cfg.trialdef.pre  * hdr.Fs);
posttrig =  round(cfg.trialdef.poststim *1000); % posttrig =  round(cfg.trialdef.post * hdr.Fs);
% disp("Pretrig")
% disp(pretrig)
% disp("Posttrig")
% disp(posttrig)


% look for the combination of a trigger "7" followed by a trigger "64"
% for each trigger except the last one
trl = [];
for j = 1:(length(value)-1)
trg1 = value(j);
trg2 = value(j+1);
if trg1==1 % && trg2==64
  trlbegin = sample(j) + pretrig;
  trlend   = sample(j) + posttrig;
  offset   = pretrig;
  newtrl   = [trlbegin trlend offset];
  trl      = [trl; newtrl];
  %disp(strcat(trlbegin,trlend));
end
end

% When calling ft_definetrial, you would specify
cfg.trialfun                = 'ft_trialfun_general';
cfg.headerfile              = confile;
cfg.datafile                = confile;
cfg.channel                 = 1:num_channels; % hdr.grad.label; 
cfg.continuous              = 'yes';
cfg.hpfilter                = 'yes';
cfg.hpfilttype              = 'firws';
cfg.hpfreq                  = 0.1; %FIXME should try something lower to appease reviewers 0.1
cfg.hpfiltdf                = 0.15;
cfg.hpfiltwintype           = 'blackman';
cfg.hpfiltdir               = 'onepass-zerophase';
cfg.dftfreq                 = 50; % removal line noise
cfg.trialdef.pre  = 0.1;
cfg.trialdef.post = 0.4;

% and you would call:
%cfg = ft_definetrial(cfg);
%data = ft_preprocessing(cfg);
%%
%% END ADDITIONS FOR EXTRA DETAIL
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13. Visually Inspect data for "bad" trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg                     = [];
cfg.channel             = {'MEG','-AG109','-AG052'};
cfg.colorgroups         = 'allblack';
cfg.viewmode            = 'vertical';
cfg.plotevents          = 'no'; 
cfg.preproc.hpfilter    = 'yes';
cfg.preproc.hpfreq      = 1;



% The following settings are useful for identifying EOG artifacts:
%   cfg.preproc.bpfilter    = 'yes'
%   cfg.preproc.bpfilttype  = 'but'
%   cfg.preproc.bpfreq      = [1 15]
%   cfg.preproc.bpfiltord   = 4
%   cfg.preproc.rectify     = 'yes'
%
% The following settings are useful for identifying muscle artifacts:
%   cfg.preproc.bpfilter    = 'yes'
%   cfg.preproc.bpfreq      = [110 140]
%   cfg.preproc.bpfiltord   =  8
%   cfg.preproc.bpfilttype  = 'but'
%   cfg.preproc.rectify     = 'yes'
%   cfg.preproc.boxcar      = 0.2



% Load the summary again so you can manually remove any bad trials
cfg                 = [];
cfg.method          = 'summary';
cfg.keepchannel     = 'yes';
deviant             = ft_rejectvisual(cfg, deviant);

cfg                 = [];
cfg.method          = 'summary';
cfg.keepchannel     = 'yes';
predeviant          = ft_rejectvisual(cfg, predeviant);


% Plot the Clean Data
cfg = [];
cfg.channel             = {'MEG','-AG109','-AG052'};
cfg.colorgroups         = 'allblack';
cfg.viewmode            = 'vertical';
cfg.preproc.hpfilter    = 'yes';
cfg.preproc.hpfreq      = 1;
cfg.plotevents          = 'no'; 
%ft_databrowser(cfg,deviant);       % Lance had error: Unrecognized field name "offset".
%ft_databrowser(cfg,predeviant);    % Lance had error: Unrecognized field name "offset".
% Save to File
save deviant deviant
save predeviant predeviant
%ft_databrowser(cfg,data_clean);     % Lance had error: Unrecognized field name "offset"





%%

%% Step 14. Simple ERF Analysis
% Some very simple ERF analysis:

cfg                 = [];
cfg.layout          = lay;
cfg.linewidth       = 2;
figure; ft_multiplotER(cfg, ft_timelockanalysis([],deviant), ...
   ft_timelockanalysis([],predeviant));
%%







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 14. Perform MRI-MEG Coreg with MEMES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load grad_trans and headshape_downsampled
disp('Loading grad_trans and headshape_downsampled');
load('grad_trans');
load('headshape_downsampled');

dir_name = folder
output_path = strcat(folder,'CHILD_MEMES_V2\')
folder_meg_research = strcat(code_folder, 'Macquarie-MEG-Research')
% Path to MRI Library for MEMES
path_to_MRI_library = strcat(folder_data,'\Child_MRI\new_HCP_library_for_MEMES\\') ;
% Path to MQ_MEG_Scripts
path_to_MQ_MEG_Scripts = strcat(folder_meg_research,'\MQ_MEG_Scripts-master\MQ_MEG_Scripts-master\') 
% Path to MEMES
path_to_MEMES = strcat(folder_meg_research,'\MEMES-master\MEMES-master\') 


cd(folder)
%cd(strcat(code_folder,'Roving_MMN_2020-master\Roving_MMN_2020-master\preprocessing_scripts_ME125_phase2_YSun\Others\child_MEMES'))
elpfile = dir('*.elp').name
hspfile = dir('*.hsp').name

% mrkfile =   '2872_CM_ME125_2018_01_12_INI.mrk'
all_mrks = dir('*.mrk')
all_mrks_table = struct2table(all_mrks); % convert the struct array to a table
all_mrks_by_size = sortrows(all_mrks_table, 'date', 'ascend'); % sort the table by 'size', as 15 min roving oddball size will be larger than 10 min resting state
all_mrks_by_size = table2struct(all_mrks_by_size);
mrkfile = all_mrks_by_size.name

% confile =   '2872_CM_ME125_2018_01_12_B1_denoise_rethm.con'
if adult == 1
    all_cons = dir('*tspca.con');
else
    all_cons = dir('*rethm_tspca.con');
    if length(all_cons) == 0
        all_cons = dir('*denoise_rethm.con');
    end
end
if length(all_cons) == 0
    if adult == 1
        all_cons = dir('*B3.con');
    elseif adult == 0
        all_cons = dir('*B1.con');
    end
end
if length(all_cons) == 0
    all_cons = dir('*.con');
end
if length(all_cons) == 0
    disp("NO CON FILE FOUND")
end
all_cons_table = struct2table(all_cons); % convert the struct array to a table
all_cons_by_size = sortrows(all_cons_table, 'bytes', 'descend'); % sort the table by 'size', as 15 min roving oddball size will be larger than 10 min resting state
all_cons_by_size = table2struct(all_cons_by_size);
confile = all_cons_by_size.name    % Gets first one somehow




% grad_trans,headshape_downsampled
% child_MEMES_v2(dir_name,elpfile,hspfile,confile,mrkfile,path_to_MRI_library);