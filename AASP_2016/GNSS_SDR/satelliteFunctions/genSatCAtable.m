function [ CATable ] = genSatCAtable( sats )
%GENSATCATABLE Generates CA look-up table of the used satellites for signal
%generation 
%   Inputs:     sats - matrix with informations regarding desired
%                       satellites
%   Outputs:    CATable - the look-up table for signal generation

% preallocation
d = length(sats);
CATable = zeros(d,1023);

%creation loop
for i = 1:d
    CATable(i,:) = generateCA(sats(i).PRN);
end %for 1:d

end

