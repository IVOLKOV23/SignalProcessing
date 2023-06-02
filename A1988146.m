%% Week 2 Workshop  - Visualising Neural Signals

%% 1. Install MATLAB & toolboxes
%
% The following MATLAB toolboxes are required for this subject.
%
% * Signal Processing Toolbox
% * Systems Identification Toolbox
% * Statistical and Machine Learning Toolbox
% * DSP System Toolbox
% * Image Processing Toolbox
% 
%% 2a. Set working directory (wd)
%
% * Download and save data for Workshop Week 2: datafile004.nev and
% datafile004.ns6
% https://drive.google.com/drive/folders/13bpRuhj-AoDjfXiaDQSv9NMclj3udo86?usp=sharing
% 
% * Download NPMK-5.5.2.0.zip from https://github.com/BlackrockMicrosystems/NPMK
% * Setup a working directory with folders for the NPMK toolbox, signal data, assigment code & reports.
%
%% 2b. A comment on commenting
% 
% For all workshops and assignments, ensure you use comments explaining each step in your code. Please include
% axis labels with units, a title, and a legend if appropriate. All figures in workshop questions are numbered
%  e.g. Figure 1a - be sure to put the figure numbering in the title of every figure. Take note
% of the scaling of each graph given in the question.
%
%% 3. Install NPMK
% installNPMK
%% 4. Signals
% In the zip folder you will find data that contain the following types of bio-signals:
%
% * Multielectrode neural recordings
% * Triggers 
%
% Set data path and load files using the code below.

% set wd
dataPath = 'C:\Users\vid15\Documents\MATLAB\Signal Processing\NPMK-master\NPMK';
% dataPath = 'C:\Users\hmeffin\Home\Teaching\2021\BMEN90035 Biosignal Processing\Workshops\Matlab\Data';
datafile = 'datafile004';

%  generate path to the NEV data file containing triggers and spike times
pathNEV = fullfile(dataPath, [datafile, '.nev']); 

%  generate path to the ns6 data file containing containing raw data
pathNS6 = fullfile(dataPath, [datafile, '.ns6']);               

% load files
NEV = openNEV(pathNEV, 'read', 'nosave');
ns6 = openNSx(pathNS6, 'read','uv');

%% 
% .nev files contain
%
% * trigger times on digital and analog channels
% * spike times according to Blackrock acquisition system (Cereplex Plus)
% * spike waveforms according to Blackrock acquisition system
% 
% .ns6 file contains
% 
% * raw voltage traces on 32 channels of electrodes
%
%% 5. Extracting data
% Extract data from the .nev, and .ns6 files.

% sampling frequency of ns6 file (Hz)
Fs6           = ns6.MetaTags.SamplingFreq; 

% one channel to read (1-32 are electrodes)
iChannel      = 4;

% voltage signal for channel iChannel (uv)
voltageSignal = ns6.Data(iChannel,:);                           

% trigger times on digital channel (sec), digital channel (sample), and analog input (samples)
TimesTriggerDigital     = NEV.Data.SerialDigitalIO.TimeStampSec;
TimesTriggerDigitalInd  = NEV.Data.SerialDigitalIO.TimeStamp;
TimesTriggerAnalogInd   = NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==129);

% spike times
spikeTime = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==iChannel))/Fs6; 

% a spike index
iSpike    = 1001;

% waveform of spike extracted by Blackrock
waveform  = double(NEV.Data.Spikes.Waveform(:,iSpike));                                 

%% 6. Visualising data: Figure 1 (4 marks)
%
% * Figure 1: Plot a subsample (every 10 samples) of the raw voltage traces for all electrodes 1-32 staggered
% vertically for visibility for the whole recording time (seconds). Use as
% staggered offset of 500 units on the y-axis between traces. Use
% different colours for different traces to further improve visibility.
% Remember by subsampling you change the sampling frequency of the data.

% all electrode channels
all_channels = 1:32;

% sub-sampled voltage signals
sub_sampled_voltageSignal = ns6.Data(all_channels,1:10:end);

% new sampling rate
new_Fs6 = Fs6/10; 

% figure 1
figure;

% x-axis (the whole recording time)
t_1 = (1:length(sub_sampled_voltageSignal))./(Fs6/10);

% plotting each channel data with an offset of 500 (1st channel at the top)
hold on
for j_1 = 1:32
    plot(t_1, sub_sampled_voltageSignal(j_1, :) - 500*(j_1 - 1));
end
hold off

% figure annotation
set(gca, 'YTick', []);
xlabel('Time (seconds)');
ylabel('Channels');
title('Figure 1 - Raw, sub-sampled voltage data from the 32 electrodes');
xlim([0 450]);

%% 7. Visualising data: Figure 2 (11 marks)
%
% Figure 2: Plot raw voltage traces for all electrodes 1-32 staggered
% vertically for visibility for 30s, 3s, and 0.3s windows (Fig 2a, 2b & 2c),
% starting at sample index 900,000. 
%
% * Figure 2a: On the plot with time window of 30s overlay the
% digital trigger times as vertical lines. Use an offset of 500 units on the y axis between
% traces and different colours for the different traces.

% raw voltage data for a 30s window
% use sample/sampling freq = time for voltageSignal calculations
raw_voltageSignal_30 = ns6.Data(all_channels, 900000:1800000-1);

% figure 2a
figure;

% x-axis (time window of 30s to 60s)
t_2a = 30 + (1:length(raw_voltageSignal_30))./(Fs6);

% plotting each channel voltage data with an offset of 500 (1st channel at the top)
hold on
for j_2a = 1:32
    plot(t_2a, raw_voltageSignal_30(j_2a, :) - 500*(j_2a - 1));
end
hold off

% digital trigger times overlay for the 30s window
xline(TimesTriggerDigital(3:8), '-');

% figure annotation
set(gca, 'YTick', []);
xlabel('Time (seconds)');
ylabel('Channels');
title('Figure 2a - Raw voltage data, 30s window, digital trigger times overlayed');
xlim([30 60]);

% * Figure 2b: Use an offset of 500 units on the y axis between
% traces and different colours for the different traces.

% raw voltage data for a 3s window
% use sample/sampling freq = time for voltageSignal calculations
raw_voltageSignal_33 = ns6.Data(all_channels, 900000:990000-1);

% figure 2b
figure;

% x-axis (time window of 30s to 33s)
t_2b = 30 + (1:length(raw_voltageSignal_33))./(Fs6);

% plotting each channel voltage data with an offset of 500 (1st channel at the top)
hold on
for j_2b = 1:32
    plot(t_2b, raw_voltageSignal_33(j_2b, :) - 500*(j_2b - 1));
end
hold off

% figure annotation
set(gca, 'YTick', []);
xlabel('Time (seconds)');
ylabel('Channels');
title('Figure 2b - Raw voltage data, 3s window');
xlim([30 33]);

% * Figure 2c: On the plot with time window of 0.3s, overlay the spike 
% times by superimposing dots in a different colour on the 
% correct voltage channel trace. Suggest you use dots as magneta ('m.')
% and 'MarkerSize',10. Use an offset of 200 units on the y axis between
% traces and black for all the traces to make the spike dots visible.
%

% raw voltage data for a 0.3s window
% use sample/sampling freq = time for voltageSignal calculations
raw_voltageSignal_303 = ns6.Data(all_channels, 900000:909000-1);

% figure 2c
figure;

% x-axis
t_2c = 30 + (1:length(raw_voltageSignal_303))./(Fs6);

hold on
for j_2c = 1:32
    
    % creating a spikeTime vector for each of the channels
    spikeTime_2c = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==j_2c))/Fs6;
    
    % plotting each channel voltage data with an offset of 200 (1st channel at the top)
    plot(t_2c, raw_voltageSignal_303(j_2c, :) - 200*(j_2c - 1), 'Color', 'k');
    
    % extracting voltage data for each channel
    voltageSignal_2c = ns6.Data(j_2c, :);
    
    for z = 1:length(spikeTime_2c)
        
        % for each spikeTime in a spikeTime vector extract only the ones
        % that lie in the specified window
        if spikeTime_2c(z) >= 30 && spikeTime_2c(z) <= 30.3
            
            % save the vector index of that spikeTime
            index = round(spikeTime_2c(z)/(1/Fs6));
            
            % superimpose spikeTimes with traces for each of the channels
            plot(spikeTime_2c(z), voltageSignal_2c(index) - 200*(j_2c - 1), '.', 'MarkerSize', 10, 'Color', 'm');
        end
        
    end
    
end
hold off

% figure annotation
set(gca, 'YTick', []);
xlabel('Time (seconds)');
ylabel('Channels');
title('Figure 2c - Raw voltage data, 0.3s window, spikeTimes superimposed');
xlim([30 30.3]);

%% 8. Interpretating data (10 marks)

% * What are the important neural signals present in these recordings
% and what are possible sources of noise?
%%
% 
% * As for the neural signals, I think one can distinguish the shapes of
% alpha, beta and maybe even theta waves, though the signal is extremely
% noisy
% * Possible sources of noise include:
% * Environmental: AC powerline interference, lighting & electronic
% equipment
% * Physiological: cardiac signals, motion artifacts, ocular signals, breathing, 
% other neural brain activity, skin/surface potentials 

%%

% * What periodicities are evident in the signal? Suggest possible interpretations: which
% periodicities might be related to neural signals and which periodicties
% might be related to noise? 

%%
% * At each of the trigger times there appears to be a spike in either
% negative or positive direction. This is certaintly related to brain
% activity - reaction to a visual stimulus
% * Before every trigger time there also appears aquick but noticable
% increase in frequency and amplitude in voltage. I am not sure if that is
% due to noise or response to a stimulus, however, I am inclined to believe
% that this is due noise, due to the erratic nature
% * There also appears sinusoidal or almost "flat" traces for some of the
% electrodes, meaning of which I am not entirelly sure. These are the only
% traces that do not seem to deep during trigger times, which leads me to
% two hypotheses. Either these are neural signals that do not respond to
% the stimuli or are sources of noise, line powerline interference due to
% the sinusoidal nature of the trace
%%

% * What events are evident in the signal (i.e. sudden changes in the signal)?  
% Suggest possible interpretations: which  events might be related to neural 
% signals and which events  might be  related to noise?  
%%
% 
% * Non-periodically there appear sudden, sharp straightlly-down vertical
% falls in some of the traces. I hypothesise that this would be due to noise
% or some sort of interference
% * At around 250 seconds there appears to be a big spike for almost all of the traces
% The source is unclear though. It is not periodic, so cannot be related to
% the brain activity stimulated by the stimuli and very sharp and quick for
% noise. It might be some other brain activity, unrelated to the stimuli
% * Innitial irregular spike (right at the beginning of the EEG recordings) 
% can be attributed to the initialisation of the set up, EEG reading and 
% switching on of the equipment
% * From time to time there are some high-amplitude spikes. It is
% impossible to determine which trace they caome from, but it appears to be
% a congregation of multiple spikes. These are probably from some sharp
% brain activity occuring during the EEG, although the voltage is pretty
% high. It cannot be breathing or cardica related, since it is not
% periodic

%% 9. Assignment Submission 
%
% BMEN90035WorkshopWeek2.m from canvas will serve as both your code and report file. Write your
% solutions in this file format and then publish this to a .pdf file. 
% Please submit the following:
%
% * The matlab code file ‘.m’ with file name ‘A1studentno.m’.
% * The published matlab report with file name ‘A1studentno.pdf’. Use the 
% inbuilt publish function in matlab to create a published version 
% of the matlab code file. 
% * A single file of any scanned mathematical derivations/solutions with file name ‘A1mstudentno.pdf’
 

