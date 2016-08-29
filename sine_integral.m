%% sine_integral
% calculates the sine integral

%% input parameters
% input: sine integral parameter

%% output parameters
% output: sine integral value

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
function output=sine_integral(input)

% allocate outputs
output=zeros(size(input));

%% approximation by numerical integration

for i=1:length(input)

    if (abs(input(i)) > eps)
        output(i)=quadgk(@sine_int,0,input(i), 'AbsTol',1.0e-10,'RelTol',0.0,'MaxIntervalCount',10000); 
    else
        output(i)=0; 
    end
    
end


%% sine integral function
function y=sine_int(x)

y=sin(x)./x;

y(isnan(y))=0;
