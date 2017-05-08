function perf = comp_perf(tgt,ref,fs)

% Cepstral distance parameters
P.cd = struct('frame' , 0.064   , ...
    'shift' , 0.016    , ...
    'window', @hanning, ...
    'order' , 24      , ...
    'timdif', 0.0     , ...
    'cmn'   , 'y');

% Frequency-weighted segmental SNR parameters
P.fwsegsnr = struct('frame'  , 0.064, ...
    'shift'  , 0.016, ...
    'window' , @hanning, ...
    'numband', 23);

% Log likelihood ratio
P.llr = struct('frame' , 0.064, ...
    'shift' , 0.016, ...
    'window', @hanning, ...
    'lpcorder', 12);


% Cepstral distance
[cd_mean, ~, timdif] = cepsdist_unsync(ref, tgt, fs, P.cd);
perf.cd = cd_mean;
if timdif < 0
    ref = ref(1 - timdif : end);
    tgt = tgt(1 : end + timdif);
else
    ref = ref(1 : end - timdif);
    tgt = tgt(1 + timdif : end);
end
% LLR
perf.llr = lpcllr(tgt, ref, fs, P.llr);
% Normalzied SRMR
perf.srmr = 10*log10(SRMR(tgt, fs, 'norm', 1));
% Frequency-weighted segemental SNR
[snr_mean] = fwsegsnr(tgt,ref,fs, P.fwsegsnr);
perf.snr = snr_mean;
%  PESQ
perf.pesq = PESQWrapperDos(ref,tgt,fs);
delete('pesq_results.txt')