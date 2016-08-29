function muest = standard_esprit_Rd(Us,M,UseJointdiag)

% STANDARD_ESPRIT_RD   R-D Standard ESPRIT algorithm
% 
% Syntax:
%   MUEST = STANDARD_ESPRIT_RD(Us,M[,UseJointdiag])
%
% Input:
%   Us - signal subspace matrix of size PROD(M) x d
%   M  - length-R vector, where M(r) isthe number of sensors in the r-th
%      mode for r=1, 2, ..., R.
%   UseJointdiag - 1 (0) if joint diagonalization of Psi-matrices should
%      (not) be performed. This parameter is optional and defaults to 1 (yes).
%      If joint diagonalization is not used, the pairing of the sources
%      across the modes will not be guaranteed.
%
% Output:
%   MUEST - matrix of size R x d containing the estimated spatial
%      frequencies for all d sources and all R modes.

R = length(M);

if nargin < 3
    UseJointdiag = 1;
end

psihat = cell(1,R);
% Solve shift invariance equations mode-by-mode
for r = 1:R
    % selection matrices for mode r
    J2 = kron(kron(speye(prod(M(1:r-1))),sparse([zeros(M(r)-1,1),eye(M(r)-1)])), speye(prod(M(r+1:R))));
    J1 = fliplr(flipud(J2));
    % shift invariance equation for mode r: solved using simple LS solution (pseudo-inverse)
    psihat{r} = pinv(J1*Us)*(J2*Us);
end

if UseJointdiag
    % Using joint diagonalization
    psihat = jointdiag_c(psihat);
    % Now the eigenvalues are on the diagonals of psihat
    for r = 1:R
        % the eigenvalues are e^(j mu_i^(r)) -> ANGLE recovers mu_i^(r)
        muest(r,:) = angle(diag(psihat{r})).';
    end
else
    % Now using joint diagonalization
    for r = 1:R
        % eigenvalues computed individually
        muest(r,:) = angle(eig(psihat{r})).';
    end
end
