%% MEG data

adults = 1
if adults~=1
    %Children
    meg_path = 'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Child_MEG\2872'
    meg_file = '2872_CM_ME125_2018_01_12_B1_denoise_rethm.con'  % your exported con-filename

    

    % %Children
    audio_ch        = 135; % Audio event from fOMRI
    every_trigger   = 145  % spikes every event (every 500ms)
    trig_chans      = 146:152 % 7-tone
    other_ch        = 192  % Unsure what this is, roughly is on/off while any event is occurring

end


if adults == 1
    % Adults
    meg_path = 'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Adult_MEG\2737\'
    % meg_file = '2737_KY_ME125_2017_10_20_B1.con';  % Resting state
    % meg_file = '2737_KY_ME125_2017_10_20_B2.con';  % Resting state
    meg_file = '2737_KY_ME125_2017_10_20_B3.con'  % Roving oddball
    % meg_file = '2737_KY_ME125_2017_10_20_B4.con';  % Resting state

    %%Probable event channels
    %Adults
    audio_ch = 166;        % Audio event from fOMRI  
    every_trigger = 193;   % spikes every event (every 500ms)
    trig_chans = 194:204   % 11 tones
    trig_chans = 194:200   % 7-tone
    other_ch   = 256  % Unsure what this is, roughly is on/off while any event is occurring
    
end


% chan_analog = 192; % Unknown information in this channel (result usually 0, occasionally 5) 
% % See meetings
% %trig_chans = 129:134 # Noise category A
% %trig_chans = 136:144 # noise category A
% %trig_chans = 145:152 # noise category B
% %trig_chans = 153:159 # noise category A
% %trig_chans = [184 185 186];
%%
dataset = fullfile(meg_path, meg_file);
%%
event = ft_read_event(dataset, 'dataformat','yokogawa_con', 'threshold', 1.6, 'chanindx', trig_chans, 'detectflank', 'up')% , ); % , 
%%
hdr = ft_read_header(dataset, 'dataformat','yokogawa_con');
%%
% search for "trigger" events
plot([event.sample], [event.value], '.')
disp(unique({event.type}))
%%
value  = [event(find(strcmp('trigger', {event.type}))).value]'
%%
sample = [event(find(strcmp('trigger', {event.type}))).sample]';
disp(value)
disp(sample)
%%
% Read output from fieldtrip
data = ft_read_data(dataset)
% % Read output from MNE
% data = ft_read_data('F:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Child_MEG\2872\2872_CM_ME125_2018_01_12_B1_denoise_rethm_raw_tsss.fif')
%%
% Examine audio bit
duration = 15;
start_secs = 25;
end_secs = start_secs + duration;
dat = ft_read_data(dataset, 'begsample', start_secs*1000, 'endsample', end_secs*1000, 'chanindx', audio_ch);
plot(dat)
%%
% Examine every trigger bit
dat = ft_read_data(dataset, 'begsample', start_secs*1000, 'endsample', end_secs*1000, 'chanindx', every_trigger);
plot(dat)
%%
% Examine a piece of data
dat = ft_read_data(dataset, 'begsample', start_secs*1000, 'endsample', end_secs*1000, 'chanindx', other_ch);
plot(dat)
%%
start_plotting = 25*1000 % 25 secs in is when data normally begins
for chan_num = trig_chans % [audio_ch,chan_analog] 
    for plot_length = [50] % [30 10] % Various lengths in secs
        disp(strcat(num2str(chan_num), " , ", num2str(plot_length)))
            % At the lowest zoom level, 
        figure;
        plot(data(chan_num,start_plotting:(start_plotting+plot_length*1000)));
    end
end
% -> % Unsure what chan_analog is, using every_trigger
% the low and high signal are most easily visible
% Halfway between high and low is a good threshold for detecting the event
%% More basic code

% From ME125 - Phase 2 - Step 1

cfg                         = [];
cfg.dataset                 = dataset;
cfg.headerfile              = dataset;
cfg.datafile                = dataset;
cfg.trialfun                = 'ft_trialfun_general';  

cfg.channel                 = hdr.grad.label; 
%cfg.channel                = 'AG*';
%cfg.channel                = '161'; % 19*'; %chan_analog
%cfg.channel                = 'all'
%cfg.channel                = 'TRIG*'

cfg.continuous              = 'yes';

%highpass filter
cfg.hpfilter                = 'yes';
cfg.hpfiltwintype           = 'blackman';
cfg.hpfilttype              = 'firws';          %'but'  'firws'
cfg.hpfiltdir               = 'onepass-zerophase';
cfg.hpfreq                  = 0.1;              %FIXME should try something lower to appease reviewers 0.1
cfg.hpfiltdf                = 0.15;             % highpass transition width
%cfg.hpfiltord               = 6;                % highpass filter order
cfg.dftfreq                 = 50;               % removal line noise

data                        = ft_preprocessing(cfg);

%lowpass filter
cfg                         = [];
cfg.lpfilter                = 'yes';
cfg.lpfilttype              = 'firws';              %'but'  'firws'
cfg.lpfiltwintype           = 'blackman';           % lowpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
cfg.lpfiltdir               = 'onepass-zerophase';  % filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
cfg.lpfreq                  = 30;
cfg.lpfiltdf                = 10;

%band stop filter
cfg.bsfilter                = 'yes';
cfg.bsfreq                  = [49.5 50.5];

data                        = ft_preprocessing(cfg,data); 
%%
% Convert a fieldtrip object to an MNE object
fiff_file  = 'E:\BigData\MEG\MRES\ME125_MMN_phase1_Yanan\Child_MEG\2872\FieldTrip data MAT to MNE data FIF.fif'
fieldtrip2fiff(fiff_file, data)
%%
%% Display waveforms using ft_databrowser
cfg = [];
cfg.blocksize = 10;
cfg = ft_databrowser(cfg, data);
%%
% Examine a piece of data
dat = ft_read_data(dataset, 'begsample', start_secs*1000, 'endsample', end_secs*1000, 'chanindx', every_trigger);
plot(dat)
%%
% Scan other channels, plot those channels with things detected
events = struct('channel', num2cell(chans))
cutoff_num_events = 100
for chan = trig_chans
    disp(chan)
    event = ft_read_event(dataset, 'dataformat','yokogawa_con', 'threshold', 1.6, 'chanindx', chan, 'detectflank', 'up')% , ); % , 
    events(chan).channel = event;
    length_structure = structfun(@(field) length(field),events(chan));
    if length_structure > cutoff_num_events;
        disp(strcat("Channel: ", num2str(chan), " events: ", num2str(length_structure)))
        dat = ft_read_data(dataset, 'begsample', start_secs*1000, 'endsample', end_secs*1000, 'chanindx', chan);        plot(dat)
        plot(dat)
    end
end
%%
% Scan known channel
event = ft_read_event(dataset,   'dataformat','yokogawa_con',  'threshold', 1.6, 'chanindx', trig_chans(1), 'detectflank', 'up')% , ); % , 
% event = ft_read_event(dataset, 'dataformat','yokogawa_con',  'threshold', 1.6, 'chanindx', trig_chans,  'detectflank', 'up' ) % No triggers detected, we are below threshold
%%
table = struct2table(event);
table(1:14,:)
%%
%Examine a dat piece - more complex, for if triggers are defined properly
dataset = fullfile(meg_path, meg_file);
hdr   = ft_read_header(dataset);

cfg =[];
cfg.dataset                 = dataset;
cfg.headerfile              = dataset;
cfg.datafile                = dataset;
cfg.trialfun                = 'ft_trialfun_general'; 


cfg.channel                 = hdr.grad.label; 
%cfg.channel                = 'AG*';
%cfg.channel                = '161'; % 19*'; %chan_analog
%cfg.channel                = 'all'
%cfg.channel                = 'TRIG*'

%highpass filter
cfg.hpfilter                = 'yes';
cfg.hpfiltwintype           = 'blackman';
cfg.hpfilttype              = 'firws';          %'but'  'firws'
cfg.hpfiltdir               = 'onepass-zerophase';
cfg.hpfreq                  = 0.1;              %FIXME should try something lower to appease reviewers 0.1
cfg.hpfiltdf                = 0.15;             % highpass transition width
%cfg.hpfiltord               = 6;                % highpass filter order
cfg.dftfreq                 = 50;               % removal line noise

data                        = ft_preprocessing(cfg,data);

%lowpass filter
cfg                         = [];
cfg.lpfilter                = 'yes';
cfg.lpfilttype              = 'firws';              %'but'  'firws'
cfg.lpfiltwintype           = 'blackman';           % lowpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
cfg.lpfiltdir               = 'onepass-zerophase';  % filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
cfg.lpfreq                  = 30;
cfg.lpfiltdf                = 10;

%band stop filter
cfg.headerfile              = dataset;
cfg.trialfun                = 'ft_trialfun_general'; 
cfg.bsfilter                = 'yes';
cfg.bsfreq                  = [49.5 50.5];
data                        = ft_preprocessing(cfg,data);



event = ft_read_event(dataset, 'threshold', 1.6, 'chanindx', trig_chans, 'detectflank', 'up'); % 166 % 'trigindx', audio_ch
event = ft_read_event(dataset, 'threshold', 1.6, 'chanindx', every_trigger, 'detectflank', 'up'); % 166 % 'trigindx', audio_ch


trl = [];

cfg.trialdef.eventtype = event.type;
% cfg.trialdef.eventtype  = 'analogtrig'; % use 'trial' for Yokogawa data
cfg.trialdef.eventvalue = event.value;
% cfg.trialdef.eventvalue = 1;
cfg.trialdef.prestim = 0.1;
cfg.trialdef.poststim = 0.4;
cfg.trigindx = every_trigger % 'chanindx' % chan_event chan_event ^ 192
cfg = ft_definetrial(cfg);

for i=1:length(event)
if strcmp(event(i).type, cfg.trialdef.eventtype)
  % it is a trigger, see whether it has the right value
  if ismember(event(i).value, cfg.trialdef.eventvalue)
    % add this to the trl definition
    begsample     = event(i).sample - cfg.trialdef.prestim*hdr.Fs;
    endsample     = event(i).sample + cfg.trialdef.poststim*hdr.Fs - 1;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trigger       = event(i).value; % remember the trigger (=condition) for each trial
    if isempty(trl)
      prevtrigger = nan;
    else
      prevtrigger   = trl(end, 4); % the condition of the previous trial
    end
    trl(end+1, :) = [round([begsample endsample offset])  trigger prevtrigger];
  end
end
end
%%
% %Filter
% 
% %Delete non-beeps
% for i = length(data.trialinfo):-1:1
%   data.trialinfo(i).latency = data.trialinfo(i).latency+12; 
%  
%    if strcmp(strtok(data.trialinfo(i).type, 'N'), 'DI') ~= 1
%       data.trialinfo(i) = []; 
%    else
%         [~, trigger]  = strtok(data.trialinfo(i).type, 'N'); % remember the trigger (=condition) for each trial
%         trigger = str2double(trigger(2));
%         data.trialinfo(i).tone = trigger;
%    end
% end
%%
function [trl, event] = mytrialfun(cfg)

% read the header information (including the sampling rate) and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset, 'threshold', 1.6, 'chanindx', trig_chans, 'detectflank', 'up');

% search for "trigger" events according to 'trigchannel' defined outside the function
value  = [event(find(strcmp(cfg.trialdef.trigchannel, {event.type}))).value]';
sample = [event(find(strcmp(cfg.trialdef.trigchannel, {event.type}))).sample]';

% creating your own trialdefinition based upon the events
trl = [];
for j = 1:length(value);
  trlbegin = sample(j) - round(cfg.trialdef.prestim  * hdr.Fs);
  trlend   = sample(j) + round(cfg.trialdef.poststim * hdr.Fs);
  offset   = -round(cfg.trialdef.prestim  * hdr.Fs);
  newtrl   = [ trlbegin trlend offset];
  trl      = [ trl ; newtrl];
end
end