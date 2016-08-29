function [ DoAs, A ] = ESPRIT(X, d)
%ESPRIT Estimation of Signal Parameters via Rotational Invariance Techniques
%   2-D DoA estimation
%   using Shift Invariance
%   WARNING: does not tolerate frequency aliasing
%
%   Usage:
%      DoAs = ESPRIT(X) : signal matrix (model order is assumed "M")
%      DoAs = ESPRIT(X, d) : model order
%      [DoAs, A] = ESPRIT(X): estimated (VMM) mixing matrix

[M, ~] = size(X);

if nargin < 2
    d = M;
end

% --------- SVD low-rank approximation ---------
[U, ~, ~] = svd(X);
Us = U(:, 1: d); % signal subspace

% --------- ESPRIT ---------

% Selection matrices
J1 = eye(M-1, M); % selection matrix first M - 1
J2 = J1(:, [M 1:M-1]); % selection matrix last M - 1

% Shift invariance equations
PSI = (pinv(J1 * Us)) * (J2 * Us); % Calculate PSI from subspace
phi = eig(PSI); % EVD of PSI
mu_hat = angle(phi); % Estimated spatial frequencies
DoAs = asin(mu_hat/pi)';

A = VMM(M,d,DoAs);
end
