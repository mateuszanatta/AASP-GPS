function [X,L_corr] = addnoise_colored(X0,SNR,PS,L_corr)

% ADDNOISE
%
% Add complex colored Kronecker-shaped Gaussian noise to noise-free (synthetic measurement) matrix/tensor 
%
% Syntax:
%   X = ADDNOISE_COLORED(X0,SNR[,PS])
% 
% Input:
%   X0 - noise-free measurement matrix/tensor (can be an n-dimensional
%      matrix of any size)
%   SNR - SNR in dB
%   PS - Signal power. Optional, defaults to one.
% 
% Output:
%   X - matrix/tensor of same size as X0, with ACGN added.
%
% Author:
%    Florian Roemer, TU Ilmenau
%    Stefanie Schwarz, TU München
% Date:
%    Nov 2011

if nargin < 3
    PS = 1;
end

% Compute noise-power
sigma = 10^(-SNR/20)/sqrt(2)*sqrt(PS);

% Draw ZMCSCG WHITE noise with variance 2sigma^2
Noise_white = (randn(size(X0))+ 0*1i*randn(size(X0)))*sigma;
% Compute correlation matrix
% rho = ones(length(M_antenna),1)*noise_corr;
% % Draw correlation factor for 
% [L_corr] = calc_l(rho,M_antenna);
% Draw ZMCSCG COLORED noise with variance 2sigma^2
Noise_color = L_corr*Noise_white;
% Add noise to noise-free measurement matrix/tensor
X = X0 + Noise_color;
