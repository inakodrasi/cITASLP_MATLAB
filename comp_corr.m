%% [R_smth, R_mean] = comp_corr(X,alp)
%  Compute online correlation matrices
%
%% OUTPUTS
%               R_smth     : online correlation matrix
%               R_mean     : batch correlation matrix
%
%% INPUTS
%               X          : STFT of microhpone signals
%               alp        : smoothing constant
%
%% FORMAT
%               X[M x ntime x nfreq]
%               R_smth{ntime}[M x M nfreq]
%               R_mean[M x M x nfreq]
%
%% AUTHORs: Ina Kodrasi

function [Ry_smth] = comp_corr(Y,alp)

[M, ntime, nfreq] = size(Y);


Ry_init = zeros(M,M,nfreq);
n_init = 5;

% Initialze it with the first 5 frames
for t_idx = 1:n_init
    for freq_idx = 1:nfreq
        cY = squeeze(Y(:,t_idx,freq_idx));
        cR = cY*cY';
        Ry_init(:,:,freq_idx) = Ry_init(:,:,freq_idx) + cR;
        
    end
end
Ry_init = Ry_init/n_init;
Ry_smth{1} = Ry_init;
Ry_smth{ntime} = [];


for t_idx = 2:ntime
    tmp = zeros(M,M,nfreq);
    for freq_idx = 1:nfreq
        
        cX = squeeze(Y(:,t_idx,freq_idx));
        cR = cX*cX';
        
        tmp(:,:,freq_idx) = cR;
    end
    
    Ry_smth{t_idx} =  (1-alp) * tmp + (alp) * Ry_smth{t_idx-1};    
end
