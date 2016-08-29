function [err,order] = rmse_Rd_mypairing(mu,muhat,errmode)

% RMSE_RD_MYPAIRING   Compute RMSE with pairing.
%
%   RMSE = RMSE_RD_MYPAIRING(MU,MUHAT) computes the RMSE between the true
% spatial frequencies MU and the estimated spatial frequencies MUHAT taking
% the pairing into account. MU and MUHAT are expected to be of size
% R x d. If R is different for MU and MUHAT, the smaller number is used.
% The number of sources d may also be different, in particular, MUHAT
% may have less columns than MU if not all of the paths were estimated.
% In this case, the "best matches" in MU are selected.
% The RMSE is by default a scalar quantity, computed by summing
% the MSE in the R-D spatial frequency plane over the d users and then
% taking the square root.
%
%  RMSE = RMSE_RD_MYPAIRING(MU,MUHAT,ERRMODE) allows to control the format
%  of the output variable. ERRMODE can be either
%      'total':   The total RMSE as a scalar [default].
%      'sources': A length-d vector of R-D RMSEs for each source.
%      'modes':    The RMSE per mode as the squareroot of the MSE summed 
%                 over the users.
%      'full':    A R x d matrix containing all the errors.
%
%  [RMSE,order] = RMSE_RD_MYPAIRING(MU,MUHAT) also outputs the ordering,
% such that MU(:,order) can be compared with MUHAT.

[R1,d] = size(mu);
[R2,dhat] = size(muhat);
R = min(R1,R2);
if nargin < 3
    errmode = 'total';
end

order = zeros(1,dhat);
assoc = 1:dhat;
available = 1:d;
for nouter = 1:dhat
    besterr = inf*ones(1,length(assoc));where = zeros(1,length(assoc));
    for ndhat = 1:length(assoc)
        errs = sum(abs(mu(1:R,available) - repmat(muhat(1:R,assoc(ndhat)),1,length(available))).^2, 1);
        [besterr(ndhat),where(ndhat)] = min(errs);
    end
    [bestbest,which] = min(besterr);
    if length(which)>1,which=which(1);end
    order(assoc(which)) = available(where((which)));
    assoc = [assoc(1:which-1), assoc(which+1:end)];
    available = [available(1:where(which)-1), available(where(which)+1:end)];
end
%order
err = abs(mu(1:R,order) - muhat(1:R,:)).^2;
switch errmode
    case 'full'
        err = sqrt(err);
    case 'modes'
        err = sqrt(sum(err.').');
    case 'sources'
        err = sqrt(sum(err));
    case 'total'
        err = sqrt(sum(sum(err)));
    otherwise
        error('Unknown value for parameter ERRMODE.');
end