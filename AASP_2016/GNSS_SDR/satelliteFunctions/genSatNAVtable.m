function [ NAVtable ] = genSatNAVtable( sats )
%GENSATNAVTABLE Generates look-up table of Navigational Data for desired 
%               satellites
%   Inputs:     sats - desired satellites
%   Outputs:    NAVtable - look-up table

% preallocation
d = length(sats);
NAVtable = ones(d,1500);

%creation loop --- implementation as Appendix B pg 168 of pdf
for jj = 1:d
    for i = 1:2:1500
        NAVtable(jj,i) = (-1)*NAVtable(jj,i);
    end %for i = 1:2:1500
end %for jj = 1:d

end

