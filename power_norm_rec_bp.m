%% power_norm_rec_bp
% calculates power normalization factor for a bandlimited rectangular pulse

%% input parameters
% B: bandwidth ( zB 1.023e6 )

%% output parameters
% pnorm: normalization factor

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
function pnorm = power_norm_rec_bp(B)

% calculating power normalization
pnorm=sqrt( quadgk(@spec_bp_rect_abs_q,-B,B) );

%% function for integration
function val=spec_bp_rect_abs_q(freq)

% code length
Lc=1023;

% code duration (1ms)
Ti=1.0e-3;

% chip time
Tc=Ti/Lc;

val= ( abs( sqrt(Tc) * sin(pi*Tc*freq)./(pi*Tc*freq) ) ).^2;