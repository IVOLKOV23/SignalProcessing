%% Week 4 Workshop - Fundamental filter concepts: Z-transform, Laplace-transform
% *Leaky Integrator Neuron*
%
% The leaky integrator neuron is defined by the differential equation, 
%
% $$C \frac{dv}{dt} = - \frac{v}{R} +  x(t) $$
%
% where $C$ is capacitance of the membrane, $V$ is the membrane
% potential of the neuron relative to resting potential, $R$
% is the membrane resistance, $\tau_m =RC$ is the membrane time constant  and
% $x(t)$ is the  input current to the neuron in question (from other
% neurons). 
%
%% 1. Simulate spike train input (2 marks)
% Simulate a spike train input to the leaky integrator neuron using a Poisson process. 
%
% $$x(t) = w_k \delta(t-t_k)$$
% 
% where $t_k$ are times according to the Poisson process and $w_k$ are
% input weights which are drawn from the Normal distribution with mean 0 nA and
% a variance of 50 nA (i.e. note 1 nA = 10^-9 A).
% By the definition of a Poisson process, the estimated probability of firing
% a spike during a short interval of duration T is $\lambda = f_r \times
% T$ where $f_r$ is the spike rate. 
% Use the values of 100Hz for the firing rate estimate, 100s for the 
% simulation time, and sample time step of 0.001 (T).
% Hints: the commands 'rand' and 'randn' can be used to generate uniformly
% and normally distributed random variables. Also see ClassEx3_1.m for code
% you can adapt.

% 
clear
fr      = 100;             % est. event occurance [Hz]
tSim    = 100;             % time of sim [s]
T       = 0.001;           % sim sample step [s]
nBins   = floor(tSim/T)+1; % sim sample points
nTrials = 1;               % number of simulations 
t       = 0:T:tSim;        % time of simulation, T defined above [s]

% Figure 1; Plot the spike train input using the following code:
% input weights
w = 50*10^-9*randn(1, nBins);

% lambda
lambda = fr*T;

% delta function
delta = rand(1, nBins) < lambda;

% spike train input
x = w.*delta;

%%
%
% plotting
figure(1); 
clf
plot(t, x, 'k.');
xlim([0 1]);
xlabel('Time [s]');
ylabel('Spike weight amplitude [A]');
title('Poisson spike train input');

%% 2. Continuous leaky integrator (5 marks)
% The following question requires mathematical derivation. Answers can be
% submitted as a scanned copy of hand written working, or using mathematical
% typesetting packages in Word or Latex (if you know how to use these).
%
% Starting with the equation for membrane potential of the leaky integrator neuron,
% with $\tau_m = 30ms$ apply the Laplace transform to show
%
% $$ v(t) = \frac{R}{\tau_m} \int_{0}^{t} dt' \exp(\frac{-t'}{\tau_{m}})x(t-t') $$
% 
% assuming that the membrane potential is initially at rest i.e. $v(0)=0$.
% Use your derviation to give the following for the leaky integrator neuron
%
% * transfer function 
% * impulse response function
% 
% in the continuous time domain (i.e. Laplace domain). 
%
% Hence show that the continuous time frequency response function is given
% by,
%
% $$\hat{H}(\Omega) = \frac{R}{1+j\Omega\tau_m }$$
%
% where $\Omega$ is the angular frequency radians per second.

%% 3. Discrete time leaky integrator (5 marks)
% 
% Convert the differential equation for the leaky integrator to a discrete time difference
% equation, using a time step  of $T$ and using the substitution $v(t)= v[n-1]$ on
% the right hand side of the equation and $x(t) = x[n]$.
% 
% Using your discrete equation, apply the z-transform to the system
% to give the following for the leaky integrator neuron 
%
% * transfer function
% * impulse response function: $$h[n] =\frac{TR}{\tau_m}(1-\frac{T}{\tau_m})^n u[n]$
% 
%
% Hence show that the discrete time frequency response function is given by:
% 
% $$H(e^{j\omega}) =\frac{TR}{\tau_m}\frac{1}{(1-(1-T/\tau_m) e^{-j\omega})} $$
% 
% where $\omega$ is the normalised angular frequency (i.e. normalised by
% the sampling frequency).
%
% Figure(2): Fig 2a.(subplot(2,1,1)) On the same graph plot the magnitude of the frequency response 
% in dB for the discrete and continuous time functions.Fig 2b.(subplot(2,1,2)) On the same graph plot the phase of the frequency response 
% in radians for the discrete and continuous time functions. 
% (Hint: you can use the command 'angle' to extract the phase of a complex
% number).  Use different colours and a continuous
% line ('-') for the continuous time signal and a dashed line ('--') for the discrete time signal. Use $R    = 10e6$ ohms and $\tau_m = 30e-3s$
% seconds. Use frequency in Hz from 0 Hz up to the Nyquist frequency. 
%

R    = 10e6;         % membrane resistance [ohms]
tau_m = 30e-3;        % membrane time constant [s]

f_s  = 1/T;                                  % sample frequency [Hz]
w    = pi*(0:floor(nBins/2))/(nBins/2);      %  *normalised* angular frequency
f_Hz = (f_s/2)*(0:floor(nBins/2))/(nBins/2); % frequency vec [Hz]

% freq response for continuous time function
H_c = R./(1 + i*(2*pi*f_Hz)*tau_m);

% magnitude in dB
Mag_c = 20*log10(abs(H_c));

% freq response for discrete time function
H_d = (T*R/tau_m).*(1./(1 - (1 - (T/tau_m)).*exp(-i*w)));

% magnitude in dB
Mag_d = 20*log10(abs(H_d));

%%
%
% plotting
figure(2);

% discrete and continuous time response (magnitude)
subplot(2,1,1);
plot(f_Hz, Mag_c, '-');
hold on
plot(f_Hz, Mag_d, '--');
xlabel('Frequency (Hz)');
ylabel('Magnitude [dB]');
title('Magnitude of a leaky integrator frequency response');
legend('Continuous time', 'Discrete time');

% phase of the frequency response
subplot(2,1,2);
plot(f_Hz, angle(H_c), '-');
hold on
plot(f_Hz, angle(H_d), '--');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase of a leaky integrator frequency response');
legend('Continuous time', 'Discrete time');

%%
%
% Compare your answers for discrete and continuous functions - which
% frequencies show agreement?
%
%%
%
% * In terms of the frequency response magnitude, both discrete and
% continuous functions show a strong agreement up until 67Hz. However, the
% strong discrepencies only occur after 150Hz
% * As for the frequency response phase, agreement between the functions is
% observed only till about 9Hz. After that there is a significant change in
% phase with discrete time function diverging, while continuous time
% function phase seems to converge to around negative 1.6 radians
%
%% 4. Simulate the membrane potential (4 marks)
%
% Calculate the membrane potential in response to the spike train input
% from Question 1 by convolving the input with the impulse response function 
% for both the continuous time and discrete time systems. 
%
% Figure 3: Compare the membrane potentials for the continuous and discrete systems
% by plotting on the same axes. Plot your results  using different colours and a continuous
% line ('-') for the continuous time signal and a dashed line ('--') for the discrete time signal.
% Plot your results between 0 and 1 second. (Hint: you can use the function
% 'conv' to convolve signals.) (Hint: be careful to scale your continuous time convolution appropriately to
% include the integration time step). 
%

% Parameters
R      = 10e6;      % membrane resistance [ohms]
tau_m   = 30e-3;     % membrane time constant [s]
t       = 0:T:tSim;  % time of simulation, T defined above [s]

% continuous impulse response
h_c = (R/tau_m)*exp(-t/tau_m);

% discrete impulse response
h_d = (T*R/tau_m)*(1 - (T/tau_m)).^(0:nBins-1);

% continuous time conv
conv_c = T.*conv(h_c, x);

% discrete time conv
conv_d = conv(h_d, x);

%%
%
% plotting
figure(3);

plot(t, conv_c(1:ceil(end/2)), '-');
hold on
plot(t, conv_d(1:ceil(end/2)), '--');
xlim([0 1]);
xlabel('Time (sec)');
ylabel('Voltage (V)');
title('Simulated membrane potential of a leaky integrator');
legend('Continuous time', 'Discrete time');

%%
%
% Compare the membrane potentials for the continuous and discrete systems
% by plotting on the same axes.
%
%%
%
% * From the Figure 3 it can be observed that continuous and discrete time
% systems produce the identical membrane potentials
%
%% 5. Autocorrelation and power spectral density (4 marks)
%
% Figure 4: Calculate and plot  the autocorrelation function (subplot(2,1,1)) and power spectral density of the
% discrete time membrane voltage signal (subplot(2,1,2)). For the ACF
% select an appropriate time scale on the x-axis to illustrate any temporal
% correlations in the membrane voltage.
% For the PSD use frequency in Hz from 0 Hz up to the Nyquist frequency.

% time scale
t_corr = -tSim:T:tSim;

% autocorrelation function
acf = xcorr(conv_d(1:ceil(end/2)), 'biased');

% power spectral density
psd = 10*log10(abs(fft(conv_d)/length(conv_d)).^2);

%%
%
% plotting
figure(4);

% ACF
subplot(2,1,1);
plot(t_corr, acf);
xlim([-0.1 0.1]);
xlabel('Time lag (seconds)');
ylabel('Autocorrelation');
title('ACF of membrane potential (mV^2/s)');

% PSD
subplot(2,1,2);
plot(f_Hz(1:end), psd(1:floor(nBins/2)+1));
xlabel('Frequency (Hz)');
ylabel('PSD [dB]');
title('PSD of membrane potential');

%%
%
% Why do correlations occur on this time scale?
%
%%
%
% * The highest peak for autocorrelation function always occurs at 0, as
% the autocorrelation function finds a correlation of a signal with itself,
% so it would be perfectly correlated by itself. Width of the correlation
% indicates the time frame where correlation is evident. Since a spike train input 
% is shifted by tk, which is defined by poisson process. 
% So that would explain why width of autocorrelation is 0.1, since this is
% a signal shift. This is why the width of the autocorrelation is from -0.1 to 0.1,
% because the whole signal is prefectly correlated with itself for that period  
% of time
%
%% 6. Introduction of firing threshold (3 mark)
%
% Now introduce a firing rate threshold of 30 mV (= 0.03 V) for the leaky integrator by assigning
% a spike every time there is a *positive* threshold crossing. i.e. assign a spike output
% variable s[n] = V_th if there is a spike and 0 otherwise. Use 
% your simulated membrane potential in Question 4 to calculate
% output spike times for the  discrete time system. 
% Figure 3: plot your results ontop of Figure 3 for the output spike train  
% between 0 and 1 second using green  circle ('go'). 
% Hint: you can find
% positive threshold crossings using statements like the following, where vc(n) is the signal at sample time n.: 

% if (vc(n) >= V_th) & (vc(n-1) < V_th) 
%


V_th  = 0.030;       % 30mV threshold voltage

% determine the points above or below threshold
for i = 1:length(t)
    if (conv_d(i) >= V_th) & (conv_d(i-1) < V_th)
        spike(i) = V_th;
    else
        spike(i) = 0;
    end
end

%%
%
% plotting
% plot only non-zero values
figure(3);
plot(t(spike>0), spike(spike>0), 'go');
legend('Continuous time', 'Discrete time', 'Spike time');
hold off

%% 7. Interpretation (2 marks)
%
% autocorrelation function
acf_2 = xcorr(spike, 'biased');

%%
%
% plotting
figure(5);
plot(t_corr, acf_2);
xlabel('Time (seconds)');
ylabel('Autocorrelation');
title('ACF of membrane potential (mV^2/s)');

%%
%
% * Is the output spike train a Poisson process?  Justify your answer by
% plotting the autocorrelation in Figure 5 and interpret.
%
%%
%
% * The spike train output is indeed a Poisson process, as the correlation
% lasts for exactly one sample (I mean that duration is 0.001, which is simulation sample step)
%
%% 8. Assignment Submission 
% BMEN90035WorkshopWeek4.m from canvas will serve as both your code and report file. Write your
% solutions in this file format and then publish this to a .pdf file. 
% Please submit the following:
%
% * The matlab code file ‘.m’ with file name ‘A3studentno.m’.
% * The published matlab report with file name ‘A3studentno.pdf’. Use the 
% inbuilt publish function in matlab to create a published version 
% of the matlab code file. 
% * A single file of any scanned mathematical derivations/solutions with file name ‘A3mstudentno.pdf’
