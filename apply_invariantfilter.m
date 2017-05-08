%% [YOUT] = apply_invariantfilter(w,Y)
%  Filter and sum in the STFT domain w'*Y
%
%% OUTPUTS
%               YOUT       : output signal STFT coefficients
%
%% INPUTS
%               w          : time time invariant filter
%               Y          : multichannel input STFT coefficients
%
%% FORMAT
%               YOUT[ntime x nfreq]
%               w[M x nfreq]
%               Y[M x ntime x nfreq]
%
%% AUTHORs: Ina Kodrasi

function [YOUT] = apply_invariantfilter(w,Y)

[~,ntime,nfreq] = size(Y);                              % Dimensions of Y
Yperm = permute(Y,[1 3 2]);                             % Switch Y to [M x nfreq x ntime]
% Repeat ntime tiles of w, i.e., w[M x nfreq*ntime] 
% Fold ntime tiles of Y, i.e., Yperm[M x nfreq*ntime]
% Apply dot product for each column by using .* and sum
YOUT = sum(conj(repmat(w,1,ntime)).*Yperm(:,:));        
YOUT = reshape(YOUT.',nfreq,ntime).';                   % Reshape YOUT into a [ntime x nfreq] matrix
                                                        % Transpose is required to align elements row-wise and not column-wise 