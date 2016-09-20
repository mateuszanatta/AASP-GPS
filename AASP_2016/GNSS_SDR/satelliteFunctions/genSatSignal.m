function [signal] = genSatSignal(currSat, settings, sig, samplesPerCode)
%receiver
% INPUTS:   currSat     - is the struct's position of de current working 
%                           satellite;
%           settings    - the settings data structure
%           sig         - is the previous satSignal array
%           samplesPerCode - length of arrays to be created
% OUTPUTS:  signal      - increase the sig for satSignal array 

%% Generation of CA lookup table ===========================================
caCode = generateCA(currSat.PRN);
%% Generation of unused Navigation lookup table ============================
navCode = ones(2,1500);
%% clocks initial settings =================================================
clk10_23 = 0;
clk1_023 = 0;
%% counters offset settings ================================================
contCA = [currSat.CodPhase:1023,1:(currSat.CodPhase-1)];

if (currSat.CodPhase > 50) % offset for Navigational Period Counter
    ii = floor(currSat.CodPhase/50);
    ij = mod(currSat.CodPhase,50);
    
else ij = currSat.CodPhase;
     ii = currSat.CodPhase;
end %if(countNavPer > 50)
contNavPer = ij;
contNavMess = ii;
%% Doppler Rate Deviation Allocation for VCOs ===============================

%% VCOs carriers Generation ================================================

currFreqL1 = settings.satL1freq+currSat.DoppErr;
%currFreqL2 = settings.satL2freq+currSat.DoppErr;
tt = (0:1/settings.samplingFreq:0.1); %generates 100 ms of signal
final = length(tt);
CosL1Carr = cos(2*pi*currFreqL1*tt);
SinL1Carr = sin(2*pi*currFreqL1*tt);
currSig = 0*SinL1Carr;

%% Signal Formation, as GNSS book Appendix B ==============================
for now = 1:final
    
end %for now = 1:final


signal = currSig + sig;
end

