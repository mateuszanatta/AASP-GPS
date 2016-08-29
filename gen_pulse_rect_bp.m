%% gen_pulse_rect_bp
% samples a bandlimited rectangular pulse shape

%% input parameters
% B: bandwidth ( zB 1.023e6 )
% rho: oversampling factor
% Lc: code length
% Ti: code duration
% Nc: number of chips (duration of pulse)

%% output parameters
% pulse: two-column vector, first column sample times, second column values

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
function pulse=gen_pulse_rect_bp(B,rho,Lc,Ti,Nc)

% interpolation buffer
Nc=Nc+1;

% sampling time
Ts=1/(2*B*rho);

% chip time
Tc=Ti/Lc;

% sample point
tau=(0.0:Ts:Nc*Tc-Ts)+Ts;
tau=[fliplr(-tau) 0.0 tau]';

% calculate pulse
pulse=sine_integral( 2*pi*B*( tau + Tc/2 ) ) - sine_integral( 2*pi*B*( tau - Tc/2 )  );

% normalize to unit power
pulse=pulse*(1 / (power_norm_rec_bp(B)*pi*sqrt(Tc)));

pulse=[tau pulse];