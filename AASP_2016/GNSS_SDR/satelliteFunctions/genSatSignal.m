function [signal] = genSatSignal(currSat, set,...
        points, CAtable, NAVtable, Ptable)
%receiver
% INPUTS:   sat         - is the current satellite;
%           options     - the settings data structure
%           CAtable     - CA look-up table for the satellite
%           NAVtable    - Navigation data look-up table for the satellite
%           Ptable      - P(y) code look-up table for the satellite
% OUTPUTS:  signal      - formed section of signal outputed from this satellite 

%% clocks/counters initial settings =======================================
% the carrier freq is about 1.5 GHz, the fastest clock; but there's no need
% to generate it and since it demands lots of memory; 

%masterCont = 1; % global counter, "the Carrier clock" 
contP10_23MHz = currSat.CodPhase*10; % P Code Clock Counter: 1 inc = 154 masterCont
contCA1_023MHz = currSat.CodPhase; % C/A Code Clock Counter: 1 inc = 10 ContP10_23Mhz
contNAV50Hz = 1; % NavCode Clock: 1 inc = 10230 ContP10_23Mhz
contNAV1500 = 1; % NavMessage Counter                             

%% Carriers generation ====================================================
% Frequency off-set
freq = (set.satL1freq + currSat.DoppErr);
tt = 0:1:(points-1);

CosL1Carr = zeros(1,points);
SinL1Carr = zeros(1,points);
%CosL1Carr = cos(2*pi*freq*tt/(set.nyquistGapgen*set.samplingFreq));
%SinL1Carr = sin(2*pi*freq*tt/(set.nyquistGapgen*set.samplingFreq));
clear tt;
% Doppler shift rate
%by now no rate, just difference

%% Signal Formation, as GNSS book Appendix B ==============================
for now = 1:points
    
end %for now = 1:points


signal = 0; %%%%%%%%%%%%%%%%%%
end % end function
