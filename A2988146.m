%% Week 3 Workshop - Frequency Representation, Sampling and Quantization

%% 1. Discrete Fourier Transform (3 marks)
% The Matlab fft command implements the Discrete Fourier Transform using the 
% Fast Fourier Transform algorithm. It does this by sampling time and frequency
% with discrete unitless indicies k (for frequency) and n (for time), as
% per the following formula
%
% $$ X[k] =\sum_{n=1}^{N}x[n] e^{-2 \pi j\frac{(k-1)(n-1)}{N}}$
%
% for $k =1,...,N$.
% If $T$ is the interval in seconds between $x[n]$ and $x[n+1]$, use the above formula,
% to  explain why, the frequency  $f$ in Hz is related to $k =  1,...,ceil(N/2)$ by
% $$ f = \frac{k-1}{NT}$.

%%
%
% * FFT is similar to DTF, but implemented in a way that would allow MATLAB to perform the operation
% * From DFT, frequency is related to k through the equation $$ f = \frac{k}{NT}$
% * However, since MATLAB cannot compute a vector from index number 0, it
% needs to start from index number 1, which changes the sum domain to 1 to +ve
% infinity, instead of from 0 (same for k)
% * To account for that change 1 must subtracted from n and k, to  preserve the actual equation. 
% This results in an equation equation $$ f = \frac{k-1}{NT}$
% * As for the domain of k, the stated domain would create a single-sided, positive FFT
% spectrum. So, only a positive half of the k samples is needed

%%
% 
% Comment on what happens when $k=ceil(N/2)+1, ..., N$ for real signals $x[n]$.
%

%%
%
% * When k values go beyond the centre-point of N/2, the FFT plotting would
% go in the negative direction till it reaches 0, and then switch direction
% again. So, a signal would mirrored in a way

%% 2. FFT of voltage trace (3 marks)
% Calculate the squared magnitute, $A^2$,  in dB, of the FFT of the voltage trace
% for channel 4 in datafile004. (Note: $D [dB] = 20 log_{10}(A) =10 log_{10} (A^2)$) 
%
% Figure 1: Visualise the magnitude of the FFT at the following resolutions, as subplots on the same graph:
%
% * subplot(2,2,1): 0-1000Hz
% * subplot(2,2,2): 0-200Hz
% * subplot(2,2,3): 0-30Hz
% * subplot(2,2,4): 0-5Hz
% 

% set wd
%cd '/Users/isabelledale/OneDrive - The University of Melbourne/BMEN90035 Biosignal Processing/Workshops'
addpath(genpath('NPMK'));
addpath(genpath('workshop_m_files_version_controlled'));

% set data path
%dataPath = '/Users/isabelledale/Documents/raw-data';
dataPath = 'C:\Users\vid15\Documents\MATLAB\Signal Processing\NPMK-master\NPMK';
datafile = 'datafile004';

%  generate path to the NEV data file containing triggers and spike times
pathNEV = fullfile(dataPath, [datafile, '.nev']); 

%  generate path to the ns6 data file containing containing raw data
pathNS6 = fullfile(dataPath, [datafile, '.ns6']);               

% load files
NEV = openNEV(pathNEV, 'read', 'nosave');
ns6 = openNSx(pathNS6, 'read', 'uv');

% sampling frequency of ns6 file (Hz)
Fs6           = ns6.MetaTags.SamplingFreq;

% voltage signal trace for channel 4 (do we divide by 4??)
voltageSignal = ns6.Data(4,:)/4;

% squared magnitude of FFT in dB for channel 4
D = 10*log10(abs(fft(voltageSignal)/length(voltageSignal)).^2);

% x-axis (frequency)
f = (0:1:length(D)-1)*Fs6/length(voltageSignal);

% plotting
figure(1);

% resolution of 0-1000Hz
subplot(2, 2, 1);
plot(f, D);
xlim([0 1000]);
xlabel('Frequency (Hz)');
ylabel('Squared magnitude (dB)');
title('FFT for 0-1000Hz');

% resolution of 0-200Hz
subplot(2, 2, 2);
plot(f, D);
xlim([0 200]);
xlabel('Frequency (Hz)');
ylabel('Squared magnitude (dB)');
title('FFT for 0-200Hz');

% resolution of 0-30Hz
subplot(2, 2, 3);
plot(f, D);
xlim([0 30]);
xlabel('Frequency (Hz)');
ylabel('Squared magnitude (dB)');
title('FFT for 0-30Hz');

% resolution of 0-5Hz
subplot(2, 2, 4);
plot(f, D);
xlim([0 5]);
xlabel('Frequency (Hz)');
ylabel('Squared magnitude (dB)');
title('FFT for 0-5Hz');

%% 3. Interpretation of spectra (9 marks)
%
% Mark the following on one of the plots that is most appropriate for the
% frequency range of each of these phenomena. (You can use the 'add sticky
% note' in Adobe Reader to mark points. Other methods are also acceptable.)
% 
% * $$\alpha$ - $$\gamma$ frequency bands, and the spectral frequency range for action
% potentials mentioned in the lectures (are we going trhough the alphabet or from the lectures)
% * the stimulus presentation frequency
% * respiration frequency, approximately 20 bpm
% * heart rate 80-180 bpm
% * power line artifact
% * mark and comment on any other peaks suggesting plausible mechanism. 
% 

%%
%
% * FFT for 0-1000Hz
% * Powerline interference at 50Hz
% * Powerline interference harmonics at 150Hz, 250Hz, 350Hz and 450Hz
% * Peak triplets at around 800Hz and 900Hz can be attributed to very fast
% ripples associated with the epileptogenic regions of hippocampus

%%
%
% * FFT for 0-200Hz
% * Powerline interference at 50Hz
% * Powerline interference harmonics at 150Hz
% * Temporal muscle contraction artifact at 60Hz
% * High frequency gamma activity at 139Hz (possibly associated with memory formation)
% * Gamma band waves are present at 30-100Hz
% * LFP are present at 1-100Hz

%%
%
% * FFT for 0-30Hz
% * Alpha waves are present at 8-13Hz
% * Beta waves are present at 13-30Hz
% * Theta waves are present at 4-8Hz
% * Action potential is present at around 10Hz

%%
%
% * FFT for 0-5Hz
% * Heart rate artifact at 1.3-3Hz
% * Respiration artifact at 0.34Hz
% * Stimulus presentation frequency at 0.22Hz
% * Delta waves are present at 0.5-4Hz
% * Heart rate artifact harmonic at 1Hz
% * Respiration artifact harmonic at 0.68Hz

%% 4. Low bandwidth transmission of the LFP for Brain Computer Interfaces (10 marks)
%
% For applications with brain computer interfaces using invasive electrode
% arrays the local field potential is often used as the signal to decode
% brain activity. One advantage is that this signal requires lower
% bandwidth, because its power spectral density lies in a comparatively low
% range. This is critical when signals must be transmitted wirelessly through
% a low bandwidth link. The required bandwidth can be reduce in two ways:
% by sampling at an appropriately low frequency, and through appropriate
% quantization of the signal. You will explore both of these methods in
% this question.
%
% *Resampling*
%
% * Resample the voltage trace signal to retain the local field potential 

% resampled signal
resampled_voltageSignal = ns6.Data(4, 1:40:end)/4;

% up to the gamma band (30-300 Hz). What is an adequate sampling frequency given these requirements?

%%
%
% For purpose of resampling to retain gamma waves adequate sampling
% frequency would be 750Hz to avoid liasing (more than twice of the highest frequency)
% so every 40 samples

%%
%
% * Figure 2: For channel 4 plot the voltage trace for the orginal and resampled
% signals between 30s and 30.3s.

% raw voltage data for a 0.3s window for channel 4
% use sample/sampling freq = time for voltageSignal calculations
raw_voltageSignal_303 = ns6.Data(4, 900000:909000-1)/4;

% resampled voltage data
resampled_voltageSignal_303 = ns6.Data(4, 900000:40:909000-1)/4;

% time domain for raw data
t_raw = 30 + (1:length(raw_voltageSignal_303))./(Fs6);

% time domain for resampled data
t_resampled = t_raw(1:40:end);

% plotting
figure(2);
plot(t_raw, raw_voltageSignal_303, '-k'); 
hold on
plot(t_resampled, resampled_voltageSignal_303, '--r', 'LineWidth', 2);
hold off
legend('Subsampled', 'Resampled');
xlabel('Time (seconds)');
ylabel('Voltage (uV)');
title('Figure 2: Raw, subsampled voltage data overlayed with the resampled data');

%%
%
% Inspecting the graph, what information has been lost by resampling at
% this rate?

%%
%
% By resampling at the frequency of 750Hz, rather than at original 30000Hz
% the signal has become smoother. On a good side, it means that noise and
% possible artifacts above 750Hz were attenuated. However, some peaks like
% at 30.03s and 30.29s were lost, which might have contained some
% valuable information about the obtained signal

%% 
% *Quantization*
%
% Using lecture notes, write a function 'lloydalgorithm.m' to implement Lylod's 
% algorithm to optimally quantize a signal. A structure for your algorithm can be found below.
% 
%%
%
% 1. Check your algorithm's accuracy by applying it to a white noise signal, 
% uniformally distributed between 0 and 1. Compare your result to the analytical solution
% given in lectures. 
%
% By applying the algorithm to a white noise R values should be $$ R(i) = \frac{D(i+1)+D(i)}{2}$
% By performing calculations according to the equation specified above,
% R(1) should equal to 0.0073, which it is indeed equal to. The equation
% holds true for all of the values in R and D vectors. Thus, it can be
% concluded that the algorithm has been written correctly

%%
%
% Specification:
%
% * Use  a 6 bit quantization. 
% * Use the change in the mean square error between iterations of the algorithm 
% to deterimine the tolerance for stoping criteria. (i.e. tol =
% abs(MSE(n)-MSE(n-1)).  Use a tol = 1.e-9;

% define the white noise signal
white_noise = rand(1, length(resampled_voltageSignal));

% * Figure 3: plot log10 of the MSE as a function of iteration.

% algorithm specifications
L_1 = 2^6;
tol_1 = 1.e-9;

% run the algorithm
[D_1, R_1, err_1, vstilde_1, MSE_1] = lloydsalgorithm(white_noise, L_1, tol_1);

% plotting
figure(3);
scatter(1:length(MSE_1), log10(MSE_1),  8, '*k');
xlabel('Iterations');
ylabel('Log10(MSE)');
title('Figure 3: MSE as a function of iterations for white noise');
grid on

% * Figure 4: plot the distributions of quantized samples of (1) the output of 
% your algorithm, and (2) the analytical solution. Use the code given below
% for plotting the distribution.

% plotting
figure(4);

% quantised signal
subplot(2, 1, 1);
histogram(vstilde_1, D_1, 'Normalization', 'pdf');
xlabel('Quantisation');
ylabel('Pdf');
title('Quantised distribution for white noise');

% input signal
subplot(2, 1, 2);
histogram(white_noise, L_1, 'Normalization', 'pdf');
xlabel('Quantisation');
ylabel('Pdf');
title('Analytical distribution for white noise');

%%
%
% 2. Run this algorithm for the voltage trace of the resampled LFP signal 
% you calculated above (i.e resampled up to the  gamma band).
% Specification:
%
% * Use  a 6 bit quantization. 
% * Use the change in the mean square error between iterations of the algorithm  
% to deterimine the tolerance for stoping criteria. (i.e. tol =
% abs(MSE(n)-MSE(n-1)).  Use a tol = 1.e-9*(max -min), where max is the
% maximum value of the signal and min is the minimum value of the signal.

% algorithm specifications
L_2 = 2^6;
tol_2 = 1.e-9*(max(resampled_voltageSignal) - min(resampled_voltageSignal));

% run the algorithm
[D_2, R_2, err_2, vstilde_2, MSE_2] = lloydsalgorithm(resampled_voltageSignal, L_2, tol_2);

% * Figure 5:  plot log10 of the MSE as a function of iteration.

% plotting
figure(5);
scatter(1:length(MSE_2), log10(MSE_2), 8, '*k');
xlabel('Iterations');
ylabel('Log10(MSE)');
title('Figure 5: MSE as a function of iterations for the resampled voltage signal');
grid on

%%
% * Figure 6: plot the  distributions of quantized samples for (1) your algorithm, 
% (2) using uniformly distrubuted bondaries.

% plotting
figure(6);

% quantised signal
subplot(2, 1, 1);
histogram(vstilde_2, D_2, 'Normalization', 'pdf');
xlabel('Quantisation');
ylabel('Pdf');
title('Quantised distribution for channel 4 voltage trace');

% input signal
subplot(2, 1, 2);
histogram(resampled_voltageSignal, L_2, 'Normalization', 'pdf');
xlabel('Quantisation');
ylabel('Pdf');
title('Analytical distribution for channel 4 voltage trace');

%%
%
%  * How does the optimal quantization differ from the one using uniformally
% distributed boundaries?
%
% Optimal quantisation for the voltage trace resembles the normal
% distribution, while for the uniformly distributed white noise, histogram
% is uniformly distributed. For the voltage trace there are also not 64
% bins present upon quntisation, indicating that not every bin is utilised
% for quantisation of the signal

%%
% For histograms in Figures 4 and 6 using the following code:

%figure(4/6)
%subplot(2,1,1)
%histogram(voltage_trace,boundaries_algorithm,'Normalization', 'pdf' )
%subplot(2,1,2)
%histogram(vs,boundaries_unifrom,'Normalization', 'pdf' )
%

%% Template code for Lloyd's algorithm
%
function [D, R, err, vstilde, MSE] = lloydsalgorithm(vs, L, tol)

% initialise levels
minvs = min(vs);
maxvs = max(vs);
R = linspace(minvs, maxvs, L) + 1e-4*randn(1, L);

% initialise boundaries
D = [minvs, zeros(1, L - 1), maxvs];

% initialise error break & iteration count
% MSE error
MSE = 0;
err = 100;

% magnitude between errors
yerr = 100;

% iteration count
iter = 0;

while yerr >= tol
    iter = iter + 1;
    err_prev = err;
    
    % calculate boundaries from the initialised levels
    for i = 1:L-1
        D(i + 1) = (R(i) + R(i + 1))/2;
    end
    
    % calculate new levels
    for j = 1:L
        
        % calculate numerator and denominator sums for R
        rnum = sum(vs(vs>=D(j) & vs<=D(j + 1)));
        rdenom = length(vs(vs>=D(j) & vs<=D(j + 1)));
        
        % specific condition for levels calculation 
        if rdenom == 0
            rnum = (D(j) + D(j + 1))/2;
            rdenom = 1;
        end
        
        % function output (new levels)
        R(j) = rnum/rdenom;
        
    end
      
    % initialise quantised signal vector
    vstilde = zeros(1, length(vs));
    
    % quantised voltage signal
    for k = 1:L
        
        % check voltage signal below D_upper and above D_low
        index = ((vs>D(k)) & (vs<=D(k+1)));
        vstilde(index) = R(k);
    
    end
    
    % errors for stopping condition
    err = sum((vs - vstilde).^2)/length(vs);
    yerr = abs(err_prev - err);
    
    % MSE for the figures
    MSE(iter) = immse(vs, vstilde);
    
end  
end 

%% 5. Assignment Submission 
%
% BMEN90035WorkshopWeek3.m from canvas will serve as both your code and report file. Write your
% solutions in this file format and then publish this to a .pdf file. 
% Please submit the following:
%
% * The matlab code file ‘.m’ with file name ‘A2studentno.m’.
% * The published matlab report with file name ‘A2studentno.pdf’. Use the 
% inbuilt publish function in matlab to create a published version 
% of the matlab code file. 
% * A single file of any scanned mathematical derivations/solutions with file name ‘A2mstudentno.pdf’
% 

