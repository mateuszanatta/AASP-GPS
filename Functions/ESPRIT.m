function [ DoAs, A ] = ESPRIT(Us, d)
%ESPRIT Estimation of Signal Parameters via Rotational Invariance Techniques
%   2-D DoA estimation
%   using Shift Invariance
%   WARNING: does not tolerate frequency aliasing
%
%   Usage:
%      DoAs = ESPRIT(Us) : signal subspace matrix (model order is assumed "M")
%      DoAs = ESPRIT(Us, d) : model order

[M, ~] = size(Us);

if nargin < 2
    d = M;
end

% --------- ESPRIT ---------

% Selection matrices
J1 = eye(M-1, M); % selection matrix first M - 1
J2 = J1(:, [M 1:M-1]); % selection matrix last M - 1

% Shift invariance equations
PSI = (pinv(J1 * Us)) * (J2 * Us); % Calculate PSI from subspace
phi = eig(PSI); % EVD of PSI
mu_hat = angle(phi); % Estimated spatial frequencies
DoAs = asin(mu_hat/pi)';

end
