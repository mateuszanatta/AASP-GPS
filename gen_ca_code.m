%% function gen_ca_code
% Generates GPS C/A code as specified in the document
%
% INTERFACE SPECIFICATION IS-GPS-200 Revision D 7 December 2004
%
% Here the code is represented at levels : -1 for 1
%                                           1 for 0
%% input parameters
% satID: satellite number (1 ... 36)
% satID 37 produces the same code as 34

%% output parameters
% G: C/A code

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
function G = gen_ca_code(satID)

%  input check
if ( length(satID) > 1) 
    error('Satellite ID must be a scalar !!!');
end

if ( (satID < 1) || (satID > 36) )
    error('Satellite ID must be between 1 and 36 !!!');
end
    
%% code parameters: two phase assignments and one code delay
cparam =[2 6 5; 3 7 6; 4 8 7; 5 9 8; 1 9 17; 2 10 18; 1 8 139; 2 9 140; 
        3 10 141; 2 3 251; 3 4 252; 5 6 254; 6 7 255; 7 8 256; 8 9 257; 9 10 258;
        1 4 469; 2 5 470; 3 6 471; 4 7 472; 5 8 473; 6 9 474; 1 3 509; 4 6 512;
        5 7 513; 6 8 514; 7 9 515; 8 10 516; 1 6 859; 2 7 860; 3 8 861; 4 9 862;
        5 10 863; 4 10 950; 1 7 947; 2 8 948];
    

%% code generation
% initial state - all ones
G1 = -1*ones(1,10);
G2 =  G1;

% select taps for G2 delay
s1 = cparam(satID, 1);
s2 = cparam(satID, 2);

% allocate
G=zeros(1023,1);

% number of samples(chips) 1023
for i=1:1023
    % Gold-code
    G(i) = G2(s1)*G2(s2)*G1(10);    % modulo-2 addition is multiplication for +/- data bits

    % generator 1 - shift reg 1
    G1 = [G1(3)*G1(10) G1(1:9)];
    
    % generator 2 - shift reg 2
    G2 = [G2(2)*G2(3)*G2(6)*G2(8)*G2(9)*G2(10) G2(1:9)];
end

%% code delay (from specification)
delay=cparam(satID, 3);
G=[G(delay+1:1023,1); G(1:delay,1)];
