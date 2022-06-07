% Creates event_meg.txt and event_meg.mat         
% Creates event_sound.txt and event_sound.mag    
SUBJECTS_TO_RUN = {'9019','9020'} % {'9002','9001'} % Comment out line if you want to run on everyone

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
type =   char(tlines(2));  % 'Child'  %'Adult'
data_root =  char(tlines(4)); % 'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Child_MEG\' 
data_root_mac = char(tlines(5));
code_folder = char(tlines(7));
code_folder_mac = char(tlines(8));
num_channels_adult =  str2double(tlines(10));
num_channels_child =  str2double(tlines(11));
audio_chan_adult =  str2double(tlines(12));
audio_chan_child =  str2double(tlines(13));
trig_chans_adult =  str2double(tlines(14));
trig_chans_child =  str2double(tlines(15));

if type == 'Adult'
    adult = 1;
    load(fullfile(data_root,'neighbours_160.mat'));
    num_channels = num_channels_adult;
    audio_chan = audio_chan_adult;
    trig_chans = trig_chans_adult:(trig_chans_adult+6)
else
    adult = 0;
    load(strcat(data_root,'neighbours_125.mat'));
    num_channels = num_channels_child;
    audio_chan = audio_chan_child;
    trig_chans = trig_chans_child:(trig_chans_child+6)    
end
addpath(genpath('E:\oDrive\OneDrive\Docs Sync\MyProjects\Code\Matlab Code'))
% If running on single person, uncomment
if ~windows_machine
    code_folder = code_folder_mac
end

folder = strcat(data_root) % folder with subject subfolders
cd(folder)


subjects = []
folder_list = workspace_find_variable_list_files(strcat(data_root,'\'));
disp(folder_list)
num_subjects = 0
for i=1:length(folder_list)
    folder = char(folder_list(i));
    regex = regexp(folder, '(\d\d\d\d)', 'match'); % Only folders containing subject codes only
    try
        subj = (regex(1))
        subjects = [subjects, subj]
    end
end
subjects = unique(subjects);

%subjects('3434') = [] % Remove them as they were recorded on adult system
if length(who('-regexp','SUBJECTS_TO_RUN')) > 0
    subjects = [SUBJECTS_TO_RUN]
end
disp(subjects)


issues_with = []
issues_with_meg = []
issues_with_sound = []
for j=1:length(subjects)
    try
        subject = subjects(j)
        folder = char(strcat(data_root,subject, '\'))
        cd(folder)
    
    
    
    
    
        % Find con file
        all_cons = dir('*.con');
        all_cons_table = struct2table(all_cons); % convert the struct array to a table
        all_cons_by_size = sortrows(all_cons_table, 'bytes', 'descend'); % sort the table by 'size', as 15 min roving oddball size will be larger than 10 min resting state
        all_cons_by_size = table2struct(all_cons_by_size);
        confile = all_cons_by_size.name    % Gets first one somehow    
        meg_path = folder; % Jordan trial  
        meg_file = confile;

    
        % Save detected MEG events as to which tones were played as a text file 
        try

            dataset = fullfile(meg_path, meg_file);
            event_meg = ft_read_event(dataset, 'dataformat','yokogawa_con', 'threshold', 1.6, 'chanindx', trig_chans, 'detectflank', 'up');% , ); % , 
            event_meg = orderfields(event_meg,{'value', 'sample', 'type', 'offset', 'duration'}); % re-order columns
            fields_to_delete = {'offset', 'duration'};
            event_meg = rmfield(event_meg,fields_to_delete);
            n=numel(event_meg);
            for k=1:n
                event_meg(k).type = str2num(event_meg(k).type) -trig_chans(1) +1;
            end
            N = numel(event_meg); 
            oldnames = {'value','sample','type'};
            newnames = {'event','time_ms','tone_number'};
            for k=1:length(oldnames)
              [event_meg(1:N).(newnames{k})] = deal(event_meg.(oldnames{k})) ;
            end
            event_meg = rmfield(event_meg,oldnames);
            % Label event
            for m = 1:length(event_meg)
                event_meg(m).event = 'trigger';
            end
            save event_meg event_meg;
            writetable(struct2table(event_meg), 'event_meg.txt');
            %%
        catch
            issues_with_meg = [issues_with_meg, subject]
        end
    
    
        %% Detect sounds
        try
            cfg                         = [];
            cfg.dataset                 = meg_file;
            cfg.path                    = meg_path;
            cfg.trialdef.prestim        = 0.1;         % pre-stimulus interval
            cfg.trialdef.poststim       = 0.4;         % post-stimulus interval
            cfg.trialdef.eventtype      = 'trigger';
            cfg.trialdef.eventvalue     = [1 7];       %LANCE UPDATED: Trigger numbers
            cfg.First_Channel           = trig_chans(1)
            cfg.Last_Channel            = trig_chans(end)
            cfg.Audio_Channel           = audio_chan
            cfg.Fs                      = 1000
        

            Subj_correctedTrig = {subject}
            %Subj_correctedTrig = {'2629'}
            L = length(Subj_correctedTrig{1});
        
            cfg.trialfun                = 'FindTriggers_AudioCh'; % 'ft_trialfun_general' %
            cfg                         = ft_definetrial(cfg);
            %alldata                     = ft_redefinetrial(cfg, data);
            event                       = cfg.event;
            trl                         = cfg.trl;
        

            %%Export sounds detected to text file. Isn't tone-specific, just shows number of repetitions
            event_sounds = event;
            N = numel(event_sounds);
            oldnames = {'type','sample','value'};
            newnames = {'event','time_ms','tone_number'}; % create oldnames' columns with newnames, then delete oldnames columns
            for k=1:length(oldnames)
              [event_sounds(1:N).(newnames{k})] = deal(event_sounds.(oldnames{k}));
            end
            delnames = {'type','value','sample','time','duration'};
            event_sounds = rmfield(event_sounds,delnames);
            for m = 1:length(event_sounds)
                event_sounds(m).event = 'sound_detected';
            end
            
            save event_sounds event_sounds;
            writetable(struct2table(event_sounds), 'event_sounds.txt');
            %%
        catch
            issues_with_sound = [issues_with_sound, subject];
        end
    catch
        issues_with = [issues_with, subject];
    end
end
disp(strcat("Issues with ", issues_with))
pwd