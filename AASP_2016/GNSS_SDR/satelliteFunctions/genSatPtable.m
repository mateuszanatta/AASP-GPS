function [ Ptable ] = genSatPtable( )
%GENSATPTABLE Generates P(y) look-up table of the used satellites for 
% signal generation 
%   Inputs:     no inputs since it's simulated with alternating 1's and
%               -1's
%   Outputs:    Ptable - the look-up table for signal generation

% preallocation

Ptable = ones(1,10230);

%creation loop
for i = 1:2:10230
    Ptable(i) = (-1)*Ptable(i);
end %for i = 1:2:10230

end