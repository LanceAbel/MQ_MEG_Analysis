% LANCE ABEL Dec 2021
% Roving oddball merge previously functional code from Hannah Rapaport into 64-bit Matlab
clearvars; close all 
random_seed = 5

run     = input('Run Number: ');                                                    % 3     (NUMBER)
% Once-per-day param set
sub     = '9002'                    % input('Subject Number: ');                    % '9003'(STRING): counts one up by 1 from 9000/9001 (trials) for each subject
inits   = 'SJ'                      % input('Subject Initials: ');                  % 'AH'  (STRING)
par_num = 1                         % input('Enter participant # child/adult: ')    % 2     (NUMBER): use 2 for 2nd adult, 3 for 3rd child

isMEGexperiment = 1;            % Are we in the MEG?
megtriglength   = 10;           % Trigger length in msecs

duration_1      = 15; %15           % Main experiment in mins  (1) skewed reps uniform tones low vol
duration_others = 10;           % Data gathering in mins   (2) skewed reps uniform tones high vol (3) uniform reps, uniform tones low vol (4) skewed reps, skewed tones low vol
% Pick experiment number and choose experiment-specfic parameters
exp_options     = [2,3,5];      % Experiments we will randomly choose from to run
if run == 1
    exp_num = 1;
    dur_exp = duration_1;
else
    % Randomised
    % exp_num = exp_options(randi(numel(exp_options), 1, 1))
    
    % Instead using pseudo-random based on D:\BHPC_Files\ME_199_Lance\New Experiment\exp_conditions.txt
    fid = fopen('D:\BHPC_Files\ME_199_Lance\New Experiment\exp_conditions.txt');
    tline = fgetl(fid)
    r = 1;
    while ischar(tline)
        if r == par_num
            exp_details = [eval(char(tline))]4
        end       
        r = r + 1;
        tline = fgetl(fid)
    end
    exp_num = exp_details(run) % 1st elem is participant #, 2nd elem is first exp_num, 3rd elem is 2nd exp_num
    fclose(fid);
    dur_exp = duration_others;
    
end
%dur_exp = 1
%exp_num = 5

clearvars r fid tline

date    = clock
year    = date(1)
month   = date(2)
day     = date(3)

experiment_folder = 'D:\BHPC_Files\ME_199_Lance\New Experiment\'
addpath(genpath(experiment_folder))

sound_folder    = strcat(experiment_folder,'\Sounds\')
stimuli_folder  = strcat(experiment_folder,'\Tone Inputs\')

data_folder     = strcat(experiment_folder,'\Exp Data\')
save_path       = strcat(data_folder,sub,'\')
try
    cd(save_path)
catch
    mkdir(data_folder,sub);
end




%% build tones
tone_low_vol = 500;
tone_high_vol = 300;

tones_low_vol = [tone_low_vol];
tones_high_vol = [tone_high_vol];
dev = [0.10 0.20 0.30 0.40 0.50 0.60]; % creates a vector with deviance from 0 to 60% from freq.
for j=1:length(dev)
    % Equivalent
    %tones_low_vol = [tones_low_vol,round(tone_low_vol*(1+dev(j)))]
    tones_low_vol = [tones_low_vol, tone_low_vol + j * 50];    % Equivalent

    % High vol linear
    % tones_high_vol = [tones_high_vol, tone_high_vol + j * 100];
end

% Non-linear high-volatility tone scheme
% tone_high_vol_additions = [200 100 50 50 100 200] % Non-linear scheme for higher-volatility tones
tones_high_vol = [300 500 600 650 700 800 1000]; 


% Set up mapping of tone IDs to frequency
dict_values_low_vol = 1:(length(dev)+1);
dict_keys_low_vol =  tones_low_vol ;
dict_values_high_vol = 1:(length(dev)+1);
dict_keys_high_vol =  tones_high_vol;
tone_freq_to_ids_low_vol = containers.Map(dict_keys_low_vol,dict_values_low_vol);
tone_freq_to_ids_low_vol(650);
tone_freq_to_ids_high_vol = containers.Map(dict_keys_high_vol,dict_values_high_vol);
tone_freq_to_ids_high_vol(650);




%% REPLACEMENT OF ABOVE CODE WITH USING FILE, uses different tones for high volatility as well
cd(sound_folder)
tone_freqs = [];
tone_ids = [];

%% Uniform tones
if exp_num == 1 % If it's the first experiment on the person, they're doing the experiment that is common across kids + adults
    input_file_name = 'min 1, max 7, skewed reps uniform tones.txt'             % Low volatility, skewed reps, uniform tones
    tone_freq_to_ids = tone_freq_to_ids_low_vol
elseif exp_num == 2
    input_file_name = 'min 1, max 7, skewed reps uniform tones high vol.txt'    % High volatility, skewed reps, uniform tones - exp 2
    tone_freq_to_ids = tone_freq_to_ids_high_vol
elseif exp_num == 3
    input_file_name = 'min 1, max 7, uniform reps uniform tones.txt'            % Low volatility, uniform reps, uniform tones - exp 3
    tone_freq_to_ids = tone_freq_to_ids_low_vol
% elseif exp_num == 4
%     input_file_name = 'min 1, max 7, uniform reps uniform tones high vol.txt'   % High volatility, uniform reps, uniform tones - no longer used
%     tone_freq_to_ids = tone_freq_to_ids_high_vol
     
%% Skewed tones
elseif exp_num == 5
    input_file_name = 'min 1, max 7, skewed reps skewed tones.txt'              % Low volatility, skewed reps, skewed tones - exp 4
    tone_freq_to_ids = tone_freq_to_ids_low_vol
end    
fid = fopen(input_file_name);


     


repeats_array = []; 
num_repeats = 0;
last_freq = 0;




while true
    thisline = fgetl(fid);
    if ischar(thisline)
        freq = str2num(thisline);
        if freq == last_freq;
            num_repeats = num_repeats + 1;
        else
            num_repeats = 1;
        end
        tone_freqs = [tone_freqs,freq]; % Keep track of all the tone frequencies in order
        tone_ids = [tone_ids,tone_freq_to_ids(freq)]; % Keep track of all the tone IDs in order
        repeats_array = [repeats_array, num_repeats]; % An array counting 1,2,3..n, 1 where there are n repeats of a tone before a deviant
        last_freq = freq;
    else
        break
    end
end
fclose(fid);


%% libraries
date_text = [num2str(year) '_' num2str(month) '_' num2str(day) ]
config_display(0, 2, [0.5, 0.5, 0.5], [0, 0, 0], 'Arial', 50, 4);
config_log(['oddball_' sub '_' inits '_' date_text '_run_' num2str(run) '_exp_' num2str(exp_num) '.log']);
config_results(['oddball_' sub '_' inits '_' date_text '_run_' num2str(run) '_exp_' num2str(exp_num) '.txt']);
config_keyboard;
%config_serial(2,19200)

config_sound;       % Configures sound
%% constant definitions
srate=44100;    % sampling rate
dur = 0.070;      % duration of each event
freq= 500;       % defines frequency of principal
dev = [0 0.10 0.20 0.30 0.40 0.50 0.60]; % creates a vector with deviance from 0 to 60% from freq.
dev = [-0.4 0 0.10 0.20 0.30 0.40 0.50 0.60 1]; % creates a vector with deviance from 0 to 60% from freq.
loud= 0.2; % loudness


%% build waves for pure tones
freqd = freq*(1+dev(:)); %    calculates new frequency
t=0:1/srate:dur;        % during 'dur' does
tones=zeros(length(dev),length(t)); %cell containing sounds
for j = 1:length(dev)
    tone = sin(2*pi*freqd(j)*t);      % creates all the wave forms to be presented and stores them in a matrix called tones
    amp = loud*tone;
    amp=wind(srate,10,amp');  % makes a ramp of 10 ms to avoid clicking
    tones(j,:) = amp'; % stores the wave forms
end


%% prepare one sound wave for 2 different output channels
% Create own wavedata 
sound_params.samprate=srate;
sound_params.t = 0:1/sound_params.samprate:0.5; %sounds 500 ms 
% Our frequencies
frequency=freq;
wave_data =sin(2*pi*frequency*sound_params.t);
% Create two rows with zeros for each of the(2) channels
no_sound_zeros=zeros(length(wave_data),1)';
% Play sound through one of the channels
channel_one_sound=[no_sound_zeros; wave_data]; % Channel 1 is silent, Channel 2 has wave_data
channel_two_sound=[wave_data; no_sound_zeros]; % Channel 2 has wave_data, Channel 1 is silent






AssertOpenGL; % Running on PTB-3? Abort otherwise.
Screen('Preference', 'SkipSyncTests', 1);
KbCheck;
WaitSecs(0.1);
GetSecs;


%% MEG: open i/o port
if isMEGexperiment
    %create IO64 interface object
    try
        p.ioObj = io64;
        % check the port status
        status = io64(p.ioObj);
    catch e
        status = 1;
        disp(['Failed to open io64: ' e.message])
    end
else
    status = 1;
end
if status == 0
    p.address = hex2dec('DFB8'); %stim2
    display(['Functioning parallel port opened at: ' num2str(p.address)])
else
    p.ioObj = [];
    p.address = [];
end

% Mke sure no triggers are active
if isMEGexperiment
    triggerstatus=1;
    while triggerstatus
        io64(p.ioObj,p.address,0);
        triggerstatus=io64(p.ioObj,p.address);
    end
end
triggerstimulus = [176,177]; % LANCE WHAT SHOULD THESE SHOW?
if isMEGexperiment
   io64(p.ioObj,p.address,triggerstimulus(1,1)-128);
end



%% Performs basic initialization of the sound driver
InitializePsychSound;
d = PsychPortAudio('GetDevices');
count = PsychPortAudio('GetOpenDeviceCount'); 



if isMEGexperiment == 1
pamaster=PsychPortAudio('Open', 14); %Panphonics
%pamaster=PsychPortAudio('Open', 1); %Panphonics
end 
if isMEGexperiment == 0
pamaster=PsychPortAudio('Open');
end



buffer_channel_one = PsychPortAudio('CreateBuffer',[], channel_one_sound);
buffer_channel_two = PsychPortAudio('CreateBuffer',[], channel_two_sound);
%PsychPortAudio('RefillBuffer', pamaster, bufferHandleC1);
PsychPortAudio('FillBuffer', pamaster, buffer_channel_one);
PsychPortAudio('FillBuffer', pamaster, buffer_channel_two);


%Lance test
sound_repetitions = 3
PsychPortAudio('Start',pamaster,sound_repetitions,0,1);
start_cogent; % initialises Matlab for running Cogent
preparestring('O',2);
%% loads sounds - writes them in buffer
for j=1:length(dev)
    loadsound([sound_folder 'tone_freq' num2str(freqd(j)) '.wav'], j);        % writes sound in the buffer with code 1, 3, 5,...
end










%% plays sounds for dur_exp
clearkeys;
%%
%%
to = time;
t = time;
c = 0;


num_stimuli_presented = 0

while(t < to + 60000*dur_exp & num_stimuli_presented < dur_exp * 60 / 0.51)
%while t < to + 60000*dur_exp % lenght of stimuli=15min*60000(ms) - this is the length in ms of the whole experiment - until the end it does
   c = c+1;              % 60000 ms is how many ms are in one minute (1s = 1000ms and 1min = 60s, so 60*1000)

   % LANCE COMMENTED OUT cogent soundplay
   % playsound(asind(c));  % plays sound with code depending on the value of the matrix "asind" in the index "i"
   
   
   freq_str = num2str(tone_freqs(c))    
   myfile = [sound_folder 'tone_freq' freq_str  '.wav']
   Fs = 44100
   [sounddata soundfreq] = audioread(myfile)
   gong = audioplayer(sounddata, Fs);
   play(gong);

   t1 = time;
   %io64(p.ioObj,p.address,asind(c));                       % Puts the pulse in
   %io64(p.ioObj,p.address,tone_freq_to_ids(freq_str));     % NEW: Puts the pulse in
   io64(p.ioObj,p.address,tone_ids(c));                     % NEW: Puts the pulse in   
   %wait(megtriglength)                                     % Pulses for a certain amount of time
   duration = 0.01
   pause(duration)                                          % Pulses for a certain amount of time
   io64(p.ioObj,p.address,0);                               % Takes the pulse out
           
    
    %%
    waituntil(t1+50);  % and waits for 100ms
    waituntil(t1+150); % outputs value 0 for acquisition computer
    waituntil(t1+500);  % and waits for 500 ms 

    t = time;
    
    %% updates logbook: if pause, tplay, and codes for freq and rep number
    %addresults(response, t1, asind(c), rsind(c));
    addresults(str2num(freq_str), t1, tone_ids(c), repeats_array(c))
    Data(c,:)=[str2num(freq_str), t1, tone_ids(c), repeats_array(c)];
    num_stimuli_presented = num_stimuli_presented + 1
end % while

%% saves Data report

%save(['oddball_Data' sub '_run' num2str(run)], 'Data')
cd(save_path)
save(['oddball_' sub '_' inits '_' date_text '_run_' num2str(run) '_exp_' num2str(exp_num)], 'Data')


PsychPortAudio('Stop', pamaster, 1);
PsychPortAudio('Close');
stop_cogent



 