function [Z] = FBA(X)
%FBA Forward-Backward Averaging
%   Doubles the amount of samples
%
%   Usage:
%      Z = FBA(X)
%

[m,n] = size(X);

R = X*X'/n;

PIM = fliplr(eye(m)); % PI-M selection matrix
PIN = fliplr(eye(n)); % PI-N selection matrix
Z = [X PIM*conj(X)*PIN]; % Forward-Backward Averaged

if ~isreal(diag(R)) || ~ishermitian(R)
    disp('WARNING: Not centrohermitian');
end

end
