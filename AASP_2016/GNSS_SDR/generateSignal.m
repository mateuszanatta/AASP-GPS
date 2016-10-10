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
if(~skip)    
    satList = input('\nSatellite Matrix = \n');
    % creating struct with satellites informations
else
    satList = SAT;
end %if(~skip)
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
end %for ii = 1:d

showSatStatus(satellites);
totalTime = now;
    disp (['-> Signal Generator started at ', datestr(totalTime)]);

    
fprintf('Processing Signal Generation: . . .\n');

%% Generation of CA look-up table =========================================
satCAtable = genSatCAtable(satellites);
fprintf('\tCA look-up table . . . . . . . . . created;\n');
%% Generation of Navigation lookup table ==================================
satNAVtable = genSatNAVtable(satellites);
fprintf('\tNavigational look-up table . . . . created;\n');
%% Generation of P_yCode - simulated as Appendix B ========================
satPtable = genSatPtable(); %not generated
fprintf('\tP(y) look-up table . . . . . . . . created;\n');
%% Generation of satellite signals ========================================
% Preallocation of satSignal array
samplesPerCode = round(settings.samplingFreq*settings.nrMSgen*...
    settings.nyquistGapgen/(settings.codeFreqBasis / settings.codeLength));
satSignal = zeros(d,samplesPerCode);

startTime = now;
    disp (['        Generation started at ', datestr(startTime)]);
if(~skip) 
    % call genSatSignal to perform the signal composition
    for ii = 1:d
        fprintf('Sat %2d ', ii)
        satSignal(ii,:) = genSatSignal(satellites(ii), settings, samplesPerCode,...
            satCAtable(ii,:), satNAVtable(ii,:), satPtable);
        fprintf(' done;\n')
    end % for ii = 1:d
else load satSignal.mat; 
end %if(~skip) 
disp(['        Generation is over (elapsed time ', ...
                                        datestr(now - startTime, 13), ')'])
fprintf('\tSatellite signals  . . . . . . . . created;\n');
%% Generation of noise for each satellite =================================


%fprintf('Noise added\n')
fprintf('\tNoise  . . . . . . . . . . . . . . added;\n');
%% Generation of multipath components =====================================
fprintf('\tMultipath components . . . . . . . added;\n');
%% Gereneration of receiver signal ========================================
% call genRcvSignal to add all signals and save in approp. file
genRcvSignal(satSignal, satellites, settings, skip);
if(skip) 
   load rcvSignal.mat; 
end %if(skip)

fprintf('\tReceiver signal  . . . . . . . . . created;\n');
disp(['-> Signal Generator is over (elapsed time ', ...
                                        datestr(now - totalTime, 13), ')'])
fprintf('\n\n');