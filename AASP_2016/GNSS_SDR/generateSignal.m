%% generateSignal.m -- starts generations of GPS data
% Autor: Alexandre Serio Buscher ==========

fprintf(['\n\n             SDR GNSS\n',...
        'Initializing Satellites Code Generation\n',...
        '----------------------------------------\n'])

%% Settings ===============================================================
% at settings struct

%% Ask inputs =============================================================
%PRN Inputs
fprintf(['\tINPUTS\n\n',...
        '\t- Remember that only 8 channel are allowed by now\n',...
        '\t- Enter a matrix where each row represents ONE satellite\n',...
        'columns:  PRN     DoA      CodePhase     DopplerError    SNR    M.Path    SIR\n',...
        '        (1 - 32)(0º-90º)(1 - 1023chips)(-10k ~ +10k)Hz (~-20dB)(l=1-LOS)(~-20dB)' ])
    
satList = input('\nSatellite Matrix = \n');
% creating struct with satellites informations
[d,~] = size(satList);
satellites = struct([]);
for ii = 1:d
    satellites(ii).PRN = satList(ii,1);
    satellites(ii).DoA = satList(ii,2);
    satellites(ii).CodPhase = satList(ii,3);
    satellites(ii).DoppErr = satList(ii,4);
    satellites(ii).SNR = satList(ii,5);
    satellites(ii).Mpath = satList(ii,6);
    satellites(ii).SIR = satList(ii,7);
end
showSatStatus(satellites);
fprintf(' Processing data: . . .\n');

%% Generation of CA look-up table =========================================
satCAtable = genSatCAtable(satellites);

%% Generation of Navigation lookup table ==================================
satNAVtable = genSatNAVtable(satellites);

%% Generation of P_yCode - simulated as Appendix B ========================
satPtable = genSatPtable(); %not generated

%% Generation of satellite signals ========================================
% Preallocation of satSignal array
samplesPerCode = round(settings.samplingFreq*settings.nrMSgen*...
    settings.nyquistGapgen/(settings.codeFreqBasis / settings.codeLength));
satSignal = zeros(d,samplesPerCode);
% call genSatSignal to perform the signal composition
for ii = 1:d
    satSignal(ii,:) = genSatSignal(satellites(ii), settings, samplesPerCode,...
        satCAtable(ii,:), satNAVtable(ii,:), satPtable);
end % for ii = 1:d

%% Generation of noise for each satellite =================================

%% Generation of multipath components =====================================

%% Gereneration of receiver signal ========================================
% call genRcvSignal to add all signals and save in approp. file

