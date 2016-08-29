%% function gen_signal
% generates a sampled gnss signal provided a chip-sequence and a pulse

%% copyright note
% Copyright (c) 2012 Manuel Stein
% Technische Universit?t M?nchen
% Department of Electrical Engineering and Information Technology
% Institute for Circuit Theory and Signal Processing
% Email : manuel.stein@tum.de
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% August 2012

%% function
function  signal=gen_signal2(tau,id,pulse,Lc,Ti,Nc)

% chip time
Tc=Ti/Lc;

% generate look up for chip index
chip_nr=repmat(mod(round(single(tau./Tc)),Lc)+1.0,1,Nc*2+1)+repmat(-Nc:1:Nc,length(tau),1);

% correct negative entires and entries bigger than Lc (circular shift)
chip_nr(chip_nr<1)=chip_nr(chip_nr<1)+Lc;
chip_nr(chip_nr>Lc)=chip_nr(chip_nr>Lc)-Lc;

% get values from code
val_code=id(int32(chip_nr));

% get sample points
samples=repmat(tau-round(tau./Tc-eps)*Tc,1,Nc*2+1)+repmat(Nc*Tc:-Tc:-Nc*Tc,length(tau),1);

% interpolate values from pulse
val_interp = interp1(pulse(:,1),pulse(:,2),samples,'cubic');

% calculate signal
% signal = sqrt(Tc)*sum(val_interp.*val_code,2);
signal=[sqrt(Tc)*sum(val_interp.*val_code,2); sqrt(Tc)*sum(val_interp.*val_code,2); sqrt(Tc)*sum(val_interp.*val_code,2); sqrt(Tc)*sum(val_interp.*val_code,2); sqrt(Tc)*sum(val_interp.*val_code,2); sqrt(Tc)*sum(val_interp.*val_code,2);...
    sqrt(Tc)*sum(val_interp.*val_code,2); sqrt(Tc)*sum(val_interp.*val_code,2); sqrt(Tc)*sum(val_interp.*val_code,2); sqrt(Tc)*sum(val_interp.*val_code,2); sqrt(Tc)*sum(val_interp.*val_code,2)];