%% Generate noiseless scenarios 

%% Scenario 1
clear all; close all; clc

% Generate RIRs
scen = 1;
fs = 16e3;
T60 = 0.36;
M = 6;
tar_idx = 1;
pos_mics = zeros(M,2);
pos_mics(:,1) = 0.08*[0:(M-1)];
str_load = strcat('Impulse_response_Acoustic_Lab_Bar-Ilan_University_(Reverberation_', ...
    num2str(T60),'0s)_8-8-8-8-8-8-8_2m_045.mat');
load(str_load)
h = resample(impulse_response,fs,48e3,100);
delay = delay_ind(h);
idx = max([1,delay-0.002*fs]);
h = h(idx:idx + T60*fs - 1,1:M);
hd = h(1:0.008*fs,:);
hr = h(0.008*fs+1:end,:);

% Generate signals
[s,fso] = audioread('sqam_mf_eng.wav');
s = resample(s,fs,fso);
x = fftfilt(h,s);
xd = fftfilt(hd,s);
xr = fftfilt(hr,s);
rms_c = rms(x(:,tar_idx));
sc_f = (1/rms_c);
x = sc_f*x;
s = sc_f*s;
xd = sc_f*xd;
xr = sc_f*xr;
perf_in = comp_perf(x(:,tar_idx),s,fs);

str_save = strcat('scenario_',num2str(scen),'_T60_',num2str(T60*1000),'ms');
save(str_save, ...
    'h', 'hd', 'hr', ...
    's', 'xd', 'xr', 'x', ...
    'pos_mics', 'perf_in')


%% Scenario 2
clear all; close all; clc

% Generate RIRs
scen = 2;
fs = 16e3;
T60 = 0.44;
M = 6;
tar_idx = 1;
pos_mics = zeros(M,2);
pos_mics(:,1) = 0.06*[0:(M-1)];
[impulse_response,fso] = audioread('Lin8Ch_503_2_RIR.wav');
h = resample(impulse_response,fs,fso,100);
delay = delay_ind(h);
idx = max([1,delay-0.002*fs]);
h = h(idx:idx + T60*fs - 1,1:M);
hd = h(1:0.008*fs,:);
hr = h(0.008*fs+1:end,:);

% Generate signals
[s,fso] = audioread('sqam_mf_eng.wav');
s = resample(s,fs,fso);
x = fftfilt(h,s);
xd = fftfilt(hd,s);
xr = fftfilt(hr,s);
rms_c = rms(x(:,tar_idx));
sc_f = (1/rms_c);
x = sc_f*x;
s = sc_f*s;
xd = sc_f*xd;
xr = sc_f*xr;
perf_in = comp_perf(x(:,tar_idx),s,fs);

str_save = strcat('scenario_',num2str(scen),'_T60_',num2str(T60*1000),'ms');
save(str_save, ...
    'h', 'hd', 'hr', ...
    's', 'xd', 'xr', 'x', ...
    'pos_mics', 'perf_in')


%% Scenario 3
clear all; close all; clc

% Generate RIRs
scen = 3;
fs = 16e3;
T60 = 0.61;
M = 6;
tar_idx = 1;
pos_mics = zeros(M,2);
pos_mics(:,1) = 0.08*[0:(M-1)];
str_load = strcat('Impulse_response_Acoustic_Lab_Bar-Ilan_University_(Reverberation_', ...
    num2str(T60),'0s)_8-8-8-8-8-8-8_2m_045.mat');
load(str_load)
h = resample(impulse_response,fs,48e3,100);
delay = delay_ind(h);
idx = max([1,delay-0.002*fs]);
h = h(idx:idx + T60*fs - 1,1:M);
hd = h(1:0.008*fs,:);
hr = h(0.008*fs+1:end,:);

% Generate signals
[s,fso] = audioread('sqam_mf_eng.wav');
s = resample(s,fs,fso);
x = fftfilt(h,s);
xd = fftfilt(hd,s);
xr = fftfilt(hr,s);
rms_c = rms(x(:,tar_idx));
sc_f = (1/rms_c);
x = sc_f*x;
s = sc_f*s;
xd = sc_f*xd;
xr = sc_f*xr;
perf_in = comp_perf(x(:,tar_idx),s,fs);

str_save = strcat('scenario_',num2str(scen),'_T60_',num2str(T60*1000),'ms');
save(str_save, ...
    'h', 'hd', 'hr', ...
    's', 'xd', 'xr', 'x', ...
    'pos_mics', 'perf_in')


%% Scenario 4
clear all; close all; clc

% Generate RIRs
scen = 4;
fs = 16e3;
T60 = 1.25;
M = 6;
tar_idx = 1;
pos_mics = zeros(M,2);
pos_mics(:,1) = 0.06*[0:(M-1)];
[impulse_response,fso] = audioread('Lin8Ch_403a_2_RIR.wav');
h = resample(impulse_response,fs,fso,100);
delay = delay_ind(h);
idx = max([1,delay-0.002*fs]);
h = h(idx:idx + T60*fs - 1,1:M);
hd = h(1:0.008*fs,:);
hr = h(0.008*fs+1:end,:);

% Generate signals
[s,fso] = audioread('sqam_mf_eng.wav');
s = resample(s,fs,fso);
x = fftfilt(h,s);
xd = fftfilt(hd,s);
xr = fftfilt(hr,s);
rms_c = rms(x(:,tar_idx));
sc_f = (1/rms_c);
x = sc_f*x;
s = sc_f*s;
xd = sc_f*xd;
xr = sc_f*xr;
perf_in = comp_perf(x(:,tar_idx),s,fs);

str_save = strcat('scenario_',num2str(scen),'_T60_',num2str(T60*1000),'ms');
save(str_save, ...
    'h', 'hd', 'hr', ...
    's', 'xd', 'xr', 'x', ...
    'pos_mics', 'perf_in')

%% Scenario 5

clear all; close all; clc
scen = 5;
fs = 16e3;
T60 = 0.73;
M = 6;
tar_idx = 1;
load('micPos_8ch.mat')
pos_mics = micPos(1:M,1:2);
[impulse_response,fso] = audioread('RIR_SimRoom3_far_AnglA.wav');
h = resample(impulse_response,fs,fso,100);
delay = delay_ind(h);
idx = max([1,delay-0.002*fs]);
h = h(idx:idx + T60*fs - 1,1:M);
hd = h(1:0.008*fs,:);
hr = h(0.008*fs+1:end,:);

% Generate signals
[s,fso] = audioread('sqam_mf_eng.wav');
s = resample(s,fs,fso);
x = fftfilt(h,s);
xd = fftfilt(hd,s);
xr = fftfilt(hr,s);
rms_c = rms(x(:,tar_idx));
sc_f = (1/rms_c);
x = sc_f*x;
s = sc_f*s;
xd = sc_f*xd;
xr = sc_f*xr;
perf_in = comp_perf(x(:,tar_idx),s,fs);

str_save = strcat('scenario_',num2str(scen),'_T60_',num2str(T60*1000),'ms');
save(str_save, ...
    'h', 'hd', 'hr', ...
    's', 'xd', 'xr', 'x', ...
    'pos_mics', 'perf_in')

