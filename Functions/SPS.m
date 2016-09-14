function [Z] = SPS(X)
%SPS Spatial Smoothing
%
%   Increases the amount of samples
%   for an observation matrix X at
%   the expense of a sensor
%
%   Usage:
%      Z = SPS(X)
%

X1 = X(1:end-1,:);
X2 = X(2:end,:);

Z = [X1 X2];

end
