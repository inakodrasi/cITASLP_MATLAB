%% Generate MDVR results for different scenarios and microphone configurations
clear all; close all;

%% Load scenario

scen = 5;
T60 = 0.73;
str_load = strcat('scenario_',num2str(scen),'_T60_',num2str(T60*1000),'ms');
load(str_load)

%% Simulation parameters
fs = 16e3;
M = 6;
% Some parameters
nFFT = 1024;                             % FFT size
nShift = 256;                            % FFT shift
winType = 'hamming';                     % FFT window
w = {'tight',winType,nFFT};
tar_idx = 1;    % Reference microphone index
c = 340;
mu = 1e-1;
Lx = length(x);
tau = 0.04;                              % Smoothing constant
alp=exp(-(nFFT-(nFFT-nShift))/fs/tau);   % Smoothing parameter
G_min = 10^(-10/10);                     % minimum gain
alp_DD = alp;
%% Generate matrix with diffuse spatial coherence

D = zeros(M,M);
for m_out = 1:M-1
    for m_in = m_out+1:M
        D(m_out,m_in) =  sqrt(  (pos_mics(m_out,1)-pos_mics(m_in,1))^2 + (pos_mics(m_out,2)-pos_mics(m_in,2))^2 ) ;
        D(m_in,m_out) = D(m_out,m_in);
    end
end
gamma = comp_diffcoh(D,nFFT,fs);

%% Compute regularized coherence matrix
nfreq = size(gamma,3);
mu = 1e-1;
for freq_idx = 2:nfreq
    gamma(:,:,freq_idx) = (1-mu)* squeeze(gamma(:,:,freq_idx)) + mu*eye(M);
end


%% Compute STFTs

for idx_m = 1:M
    X(idx_m,:,:) = dgtreal(x(:,idx_m),w,nShift,nFFT).';  
    XD(idx_m,:,:) = dgtreal(xd(:,idx_m),w,nShift,nFFT).';  
    XR(idx_m,:,:) = dgtreal(xr(:,idx_m),w,nShift,nFFT).';  
end
ntime = size(X,2);

%% Compute steering vector based on hd
d = zeros(M,nfreq);
for ind_m = 1:M
    d(ind_m,:) = freqz(hd(:,ind_m),1,nFFT/2+1,fs);
end
d_target = d(tar_idx,:);
d = d./repmat(d_target,M,1);


%% Compute MVDR filter coefficients

w_mvdr = zeros(M,nfreq);
fc = zeros(nfreq,1);
for freq_idx = 2:nfreq
    c_d = squeeze(d(:,freq_idx));
    c_coh = squeeze(gamma(:,:,freq_idx));
    fc(freq_idx) = real(( (c_d')/(c_coh)*c_d ));
    w_mvdr(:,freq_idx) = ( c_coh\c_d ) / fc(freq_idx);
end

%% Apply MVDR coefficients

[XOUT_MVDR] = apply_invariantfilter(w_mvdr,X);
% Compute performance of the MVDR
xout_mvdr = idgtreal(XOUT_MVDR.',w,nShift,nFFT,Lx);
perf_mvdr = comp_perf(xout_mvdr,s,fs);


%% Compute ML-based postfilter
[Rx] = comp_corr(X,alp);

% Estimate PSDs

phir_ml = zeros(ntime,nfreq);
phis_ml = zeros(ntime,nfreq);
phir_ml_out = zeros(ntime,nfreq);

I = eye(M);
for freq_idx = 2:nfreq
    c_d = squeeze(d(:,freq_idx));
    c_coh = squeeze(gamma(:,:,freq_idx));
    c_fc = fc(freq_idx);
    c_wmvdr = w_mvdr(:,freq_idx);
    
    for t_idx = 1:ntime
        cRx = squeeze(Rx{t_idx}(:,:,freq_idx));
        cphir = (1/(M-1)) * real( trace ( ( I - c_d*c_wmvdr' ) * cRx / c_coh ) );
        cphir = max(cphir,0);
        cphir_out = cphir/c_fc;
        
        cRxd = cRx - cphir * c_coh;      
        
        cphis = real( c_wmvdr' * ( cRxd ) * c_wmvdr );
        cphis = max(cphis,0);
        
        % Store estimated PSDs
        phir_ml(t_idx,freq_idx) = cphir;
        phir_ml_out(t_idx,freq_idx) = cphir_out;
        phis_ml(t_idx,freq_idx) = cphis;
        
    end
    
end

% Compute postfilter
rho_ml = phis_ml ./ phir_ml_out;
G_ml = rho_ml ./ (1 + rho_ml);
G_ml = max(G_ml,G_min);

%% Apply ML-based postfilter

% Apply postfilter
XOUT_ML = XOUT_MVDR .* G_ml;
xout_ml = idgtreal(XOUT_ML.',w,nShift,nFFT,Lx);

% Compute performance
perf_ml = comp_perf(xout_ml,s,fs);

%% Compute ML postfilter using DD approach

rho_ml_DD =  est_rho_DD(XOUT_MVDR,phir_ml_out,alp_DD,G_min);
% Postfilter
G_ml_DD = rho_ml_DD ./ (1 + rho_ml_DD);
XOUT_ML_DD = XOUT_MVDR .* G_ml_DD;
xout_ml_DD = idgtreal(XOUT_ML_DD.',w,nShift,nFFT,Lx);
% Performance
perf_ml_DD = comp_perf(xout_ml_DD,s,fs);

%% Compute EVD-based postfilters

phir_tr = zeros(ntime,nfreq);
phir_eig2 = zeros(ntime,nfreq);
phir_out_tr = zeros(ntime,nfreq);
phir_out_eig2 = zeros(ntime,nfreq);

for freq_idx = 2:nfreq
    c_d = squeeze(d(:,freq_idx));
    c_coh = squeeze(gamma(:,:,freq_idx));
    c_fc = fc(freq_idx);
    
    for t_idx = 1:ntime 
        
        cRx = squeeze(Rx{t_idx}(:,:,freq_idx));
        pwcRx = (c_coh \ cRx);
        [U,V] = eig(pwcRx);
        eigvals = diag(V);
        eigvals = sort(real(eigvals),'descend');
        eigvals(eigvals<0) = 0;
        cphir_tr = mean(eigvals(2:end)) ;
        cphir_eig2 = eigvals(2) ;
        cphir_tr_out = cphir_tr/c_fc;
        cphir_eig2_out = cphir_eig2/c_fc;
        
        phir_tr(t_idx,freq_idx) = cphir_tr;        
        phir_eig2(t_idx,freq_idx) = cphir_eig2;
        phir_out_tr(t_idx,freq_idx) = cphir_tr_out;        
        phir_out_eig2(t_idx,freq_idx) = cphir_eig2_out;
        
    end
end
        
% Estimate rho using the DD approach
rho_tr =  est_rho_DD(XOUT_MVDR,phir_out_tr,alp_DD,G_min);
rho_eig2 =  est_rho_DD(XOUT_MVDR,phir_out_eig2,alp_DD,G_min);
% Postfilter
G_tr = rho_tr ./ (1 + rho_tr);
G_tr = max(G_tr,G_min);
G_eig2 = rho_eig2 ./ (1 + rho_eig2);
G_eig2 = max(G_eig2,G_min);


%% Apply EVD-based postfilters
XOUT_TR = XOUT_MVDR .* G_tr;
xout_tr = idgtreal(XOUT_TR.',w,nShift,nFFT,Lx);

XOUT_EIG2 = XOUT_MVDR .* G_eig2;
xout_eig2 = idgtreal(XOUT_EIG2.',w,nShift,nFFT,Lx);

% Compute performance
perf_tr = comp_perf(xout_tr,s,fs);
perf_eig2 = comp_perf(xout_eig2,s,fs);


% %% Apply only a single-channel Wiener filter
% 
% % Estimate rho
% X_TAR = squeeze(X(tar_idx,:,:));
% 
% rho_tr_sc =  est_rho_DD(X_TAR,phir_tr,alp_DD,G_min);
% rho_eig2_sc =  est_rho_DD(X_TAR,phir_eig2,alp_DD,G_min);
% rho_ml_sc = phis_ml ./ phir_ml;
% 
% % Postfilters
% G_tr_sc = rho_tr_sc ./ (1 + rho_tr_sc);
% G_tr_sc = max(G_tr_sc,G_min);
% G_eig2_sc = rho_eig2_sc ./ (1 + rho_eig2_sc);
% G_eig2_sc = max(G_eig2_sc,G_min);
% G_ml_sc = rho_ml_sc ./ (1 + rho_ml_sc);
% G_ml_sc = max(G_ml_sc,G_min);
% 
% % Apply postfilters
% XOUT_TR_SC = X_TAR .* G_tr_sc;
% xout_tr_sc = idgtreal(XOUT_TR_SC.',w,nShift,nFFT,Lx);
% XOUT_EIG2_SC = X_TAR .* G_eig2_sc;
% xout_eig2_sc = idgtreal(XOUT_EIG2_SC.',w,nShift,nFFT,Lx);
% XOUT_ML_SC = X_TAR .* G_ml_sc;
% xout_ml_sc = idgtreal(XOUT_ML_SC.',w,nShift,nFFT,Lx);
% 
% perf_tr_sc = comp_perf(xout_tr_sc,s,fs);
% perf_eig2_sc = comp_perf(xout_eig2_sc,s,fs);
% perf_ml_sc = comp_perf(xout_ml_sc,s,fs);

%% Save results
str_save = strcat('results_',num2str(scen),'_T60_',num2str(T60*1000),'ms_M',num2str(M));
save(str_save, ...
    'xout_mvdr','perf_mvdr', 'xout_ml', 'perf_ml', 'xout_eig2', 'perf_eig2', 'xout_tr', 'perf_tr','xout_ml_DD', 'perf_ml_DD')

%% Generate plots

clear all; close all
scen = 1:5;
T60 = [0.36 0.44 0.61 1.25 0.73];
M = [2 4 6];

for idx_scen = 1:length(scen)
    cscen = scen(idx_scen);
    cT60 = T60(idx_scen);
    % Load scenario
    str_load = strcat('scenario_',num2str(cscen),'_T60_',num2str(cT60*1000),'ms');
    load(str_load)

    for idx_M = 1:length(M)
        cM = M(idx_M);
        % Load results
        str_load = strcat('results_',num2str(cscen),'_T60_',num2str(cT60*1000),'ms','_M',num2str(cM));
        load(str_load);
        
        % MVDR
        dpesq_mvdr(idx_M,idx_scen) = perf_mvdr.pesq - perf_in.pesq;
        dsnr_mvdr(idx_M,idx_scen) = perf_mvdr.snr - perf_in.snr;
        dcd_mvdr(idx_M,idx_scen) = perf_mvdr.cd - perf_in.cd;
        dsrmr_mvdr(idx_M,idx_scen) = perf_mvdr.srmr - perf_in.srmr;

        % ML
        dpesq_ml(idx_M,idx_scen) = perf_ml.pesq - perf_in.pesq;
        dsnr_ml(idx_M,idx_scen) = perf_ml.snr - perf_in.snr;
        dcd_ml(idx_M,idx_scen) = perf_ml.cd - perf_in.cd;
        dsrmr_ml(idx_M,idx_scen) = perf_ml.srmr - perf_in.srmr;
        
        % ML - DD
        dpesq_ml_DD(idx_M,idx_scen) = perf_ml_DD.pesq - perf_in.pesq;
        dsnr_ml_DD(idx_M,idx_scen) = perf_ml_DD.snr - perf_in.snr;
        dcd_ml_DD(idx_M,idx_scen) = perf_ml_DD.cd - perf_in.cd;
        dsrmr_ml_DD(idx_M,idx_scen) = perf_ml_DD.srmr - perf_in.srmr;

        % eig2
        dpesq_eig2(idx_M,idx_scen) = perf_eig2.pesq - perf_in.pesq;
        dsnr_eig2(idx_M,idx_scen) = perf_eig2.snr - perf_in.snr;
        dcd_eig2(idx_M,idx_scen) = perf_eig2.cd - perf_in.cd;
        dsrmr_eig2(idx_M,idx_scen) = perf_eig2.srmr - perf_in.srmr;
        
        % trace
        dpesq_tr(idx_M,idx_scen) = perf_tr.pesq - perf_in.pesq;
        dsnr_tr(idx_M,idx_scen) = perf_tr.snr - perf_in.snr;
        dcd_tr(idx_M,idx_scen) = perf_tr.cd - perf_in.cd;
        dsrmr_tr(idx_M,idx_scen) = perf_tr.srmr - perf_in.srmr;

    end
end

close all

idx_scen = 5;
figure
subplot(3,1,1)
bar(M,[dpesq_mvdr(:,idx_scen) dpesq_ml(:,idx_scen) dpesq_ml_DD(:,idx_scen) dpesq_tr(:,idx_scen) dpesq_eig2(:,idx_scen)])
grid on
legend('MVDR', 'ML', 'ML-DD', '\lambda_1', '\lambda_2','Orientation','Horizontal')
ylabel('\Delta PESQ')
subplot(3,1,2)
bar(M,[dsnr_mvdr(:,idx_scen) dsnr_ml(:,idx_scen) dsnr_ml_DD(:,idx_scen) dsnr_tr(:,idx_scen) dsnr_eig2(:,idx_scen)])
grid on
ylabel('\Delta fwSSNR [dB]')
grid on
subplot(3,1,3)
bar(M,[dcd_mvdr(:,idx_scen) dcd_ml(:,idx_scen) dcd_ml_DD(:,idx_scen) dcd_tr(:,idx_scen) dcd_eig2(:,idx_scen)])
grid on
ylabel('\Delta CD [dB]')
grid on
