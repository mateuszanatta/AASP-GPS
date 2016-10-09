function [signal] = genSatSignal(currSat, options,...
        points, CAtable, NAVtable, Ptable)
%receiver
% INPUTS:   sat         - is the current satellite;
%           options     - the settings data structure
%           CAtable     - CA look-up table for the satellite
%           NAVtable    - Navigation data look-up table for the satellite
%           Ptable      - P(y) code look-up table for the satellite
% OUTPUTS:  signal      - formed section of signal outputed from this satellite 

%% clocks initial settings ================================================
% the carrier freq is about 1.5 GHz, the fastest clock; but there's no need
% to generate it and since it demands lots of memory; 

%masterCont = 1; % global counter, "the Carrier clock" 
contP10_23MHz = currSat.CodPhase*10; % P Code Clock Counter: 1 inc = 154 masterCont
contCA1_023MHz = currSat.CodPhase; % C/A Code Clock Counter: 1 inc = 10 ContP10_23Mhz
contNAV50Hz = 1; % NavCode Clock: 1 inc = 10230 ContP10_23Mhz
contNAV1500 = 1; % NavMessage Counter

time = (0:1/(settings.samplingFreq*totalPoints*nyquistGap*numMS):numMS);                             
final = length(1,time);
%% counters offset settings ===============================================

%% Doppler Rate Deviation Allocation for VCOs ===============================
%by now no rate, just difference
CosL1Carr = zeros(2*pi*currFreqL1*tt);
SinL1Carr = sin(2*pi*currFreqL1*tt); 
for kk = 1:final
    freq = (settings.satL1freq+currSat.DoppErr);
    CosL1Carr = cos(2*pi*currFreqL1*tt);
    SinL1Carr = sin(2*pi*currFreqL1*tt); 
end %for kk

%% VCOs carriers Generation ================================================
%{
currFreqL1 = settings.satL1freq+currSat.DoppErr;
%currFreqL2 = settings.satL2freq+currSat.DoppErr;
tt = (0:1/settings.samplingFreq:0.1); %generates 100 ms of signal
final = length(tt);
CosL1Carr = cos(2*pi*currFreqL1*tt);
SinL1Carr = sin(2*pi*currFreqL1*tt);
currSig = 0*SinL1Carr;
%}
%% Signal Formation, as GNSS book Appendix B ==============================
for now = 1:final
    
end %for now = 1:final


signal = 0; %%%%%%%%%%%%%%%%%%
end % end function
%{
NAVtable1 = ones(2,1500);
xor = 1;
for jj = 1:1500
    xor = xor*(-1);
    NAVtable1(2,:) = xor;
end

NAVtable2 = ones(2,1500);
xor = -1;
for jj = 1:1500
    xor = xor*(-1);
    NAVtable2(2,:) = xor;
end

NAVtable3 = ones(2,1500);
xor = 1;
mor = -1;
for jj = 1:1500
    if(iseven(jj))
        xor = xor*(-1);
        NAVtable3(2,:) = xor;
    else
        mor = mor*(-1);
        NAVtable3(2,:) = mor;
    end
end
%}