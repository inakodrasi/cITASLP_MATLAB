
%% [YOUT_G,des_pow,G] = est_Phis_DD(Y_tar,int_pow,G_min,alp)
%  Estimate desired PSD using the decision directed approach and a Wiener
%  gain function
%
%% OUTPUTS        
%               YOUT_G     : output signal STFT coefficients
%               des_pow    : desired signal PSD
%               G          : Wiener gain function
%
%% INPUTS        
%               Y_tar      : target microphone signal
%               int_pow    : interference PSD (reverb+noise)
%               G_min      : minimum gain
%               alp        : smoothing parameter
%                           
%% FORMAT        
%               YOUT_G[ntime x nfreq]
%               des_pow[ntime x nfreq]
%               G[ntime x nfreq]
%               Y_tar[ntime x nfreq]
%               int_pow[ntime x nfreq]
%% AUTHORs: Ina Kodrasi


function [rho] = est_rho_DD(Y_tar,int_pow,alp,G_min)

SIR_min = 0;

ntime = size(Y_tar,1);
nfreq = size(Y_tar,2);
Y_per = Y_tar.*conj(Y_tar);
              
rho = zeros(ntime,nfreq);                        

for t_idx = 1:ntime
    % compute a posteriori SIR
    sirPost = Y_per(t_idx,:) ./ int_pow(t_idx,:);
    % compute a priori SIR
    if t_idx == 1
        sirPrio = max(0, sirPost - 1);
        %sirPrio = max(0, Y_per(t_idx,:) ./ pow_init - 1);
    else
        sirPrio = alp .* abs(XOUT).^2 ./ revPowLast + ...
            (1 - alp) .* max(0, sirPost - 1);
    end
    revPowLast = int_pow(t_idx,:);
    % threshold SIR
    sirPrio = max(sirPrio, SIR_min);
    % compute and threshold the Wiener gain
    gain = sirPrio ./ (1 + sirPrio) ;
    gain = max(gain,G_min);    
    XOUT = gain.* Y_tar(t_idx,: );
    rho(t_idx,:) = sirPrio;
end


