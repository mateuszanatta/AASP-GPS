function usedSats = initSats()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Satellites =============================================================
% Use the pre-positioned "%{" and "%}" to disable one satellite
% change the Doa, and phase values
%%{
usedSats(1).PRN = 1; % PRN/CA Code pattern. No need for change
usedSats(1).DoA = 46; % Degrees. Caution with elevation mask
usedSats(1).phase = 80; 
usedSats(1).CAcode = generateCAcode(1);
%}
%{
usedSats(2).PRN = 2; % PRN/CA Code pattern. No need for change
usedSats(2).DoA = 27; % Degrees. Caution with elevation mask
usedSats(2).phase = 512; 
usedSats(2).CAcode = generateCAcode(2);
%}
%%{
usedSats(3).PRN = 3; % PRN/CA Code pattern. No need for change
usedSats(3).DoA = 85; % Degrees. Caution with elevation mask
usedSats(3).phase = 825;  
usedSats(3).CAcode = generateCAcode(3);
%}
%%{
usedSats(4).PRN = 4; % PRN/CA Code pattern. No need for change
usedSats(4).DoA = 32; % Degrees. Caution with elevation mask
usedSats(4).phase = 330; 
usedSats(4).CAcode = generateCAcode(4);
%}
%%{
usedSats(9).PRN = 9; % PRN/CA Code pattern. No need for change
usedSats(9).DoA = 9; % Degrees. Caution with elevation mask
usedSats(9).phase = 182;  
usedSats(9).CAcode = generateCAcode(9);
%}
%%{
usedSats(12).PRN = 12; % PRN/CA Code pattern. No need for change
usedSats(12).DoA = 11; % Degrees. Caution with elevation mask
usedSats(12).phase = 904;  
usedSats(12).CAcode = generateCAcode(12);
%}
%%{
usedSats(16).PRN = 16; % PRN/CA Code pattern. No need for change
usedSats(16).DoA = 46; % Degrees. Caution with elevation mask
usedSats(16).phase = 80;  
usedSats(16).CAcode = generateCAcode(16);
%}
%%{
usedSats(22).PRN = 22; % PRN/CA Code pattern. No need for change
usedSats(22).DoA = 24; % Degrees. Caution with elevation mask
usedSats(22).phase = 200; % s. Delays between 65 to 83 ms 
usedSats(22).CAcode = generateCAcode(22);
%}
%%{
usedSats(31).PRN = 31; % PRN/CA Code pattern. No need for change
usedSats(31).DoA = 77; % Degrees. Caution with elevation mask
usedSats(31).phase = 239;  
usedSats(31).CAcode = generateCAcode(31);
%}


end

