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
else
    satList = SAT;
end %if(~skip)

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
end %for ii = 1:d

showSatStatus(satellites);
totalTime = now;
    disp (['-> Signal Generator started at ', datestr(totalTime)]);
  
fprintf('Processing Signal Generation: . . .\n');

%% Generation of CA look-up table =========================================
if(~skip)
    satCAtable = genSatCAtable(satellites);
    fprintf('\t-CA look-up table . . . . . . . . . created;\n');
end
%% Generation of Navigation lookup table ==================================
if(~skip)
    satNAVtable = genSatNAVtable(satellites);
    fprintf('\t-Navigational look-up table . . . . created;\n');
end
%% Generation of P_yCode - simulated as Appendix B ========================
if(~skip)
    satPtable = genSatPtable(); %not generated
    fprintf('\t-P(y) look-up table . . . . . . . . created;\n');
end
%% Generation of satellite signals ========================================
% Preallocation of satSignal array
samplesPerCode = round(settings.samplingFreq*settings.nrMSgen*...
    settings.nyquistGapgen/(settings.codeFreqBasis / settings.codeLength));

if(~skip)
satSignal = zeros(d,samplesPerCode);
startTime = now;
    disp (['        Generation started at ', datestr(startTime)]);
 
    fprintf('\t\t("." = 200us; "|" = 1 ms; "§" = 1 comp. C/A code)\n')
    % call genSatSignal to perform the signal composition
    for ii = 1:d
        fprintf('Sat %2d ', ii)
        satSignal(ii,:) = genSatSignal(satellites(ii), settings, samplesPerCode,...
            satCAtable(ii,:), satNAVtable(ii,:), satPtable);
        fprintf(' done;\n')
    end % for ii = 1:d
    save(cat(2,settings.path,'satSignal.mat'),'satSignal');
    disp(['        Generation is over (elapsed time ', ...
                                        datestr(now - startTime, 13), ')'])
    fprintf('\t-Satellite signals  . . . . . . . . created;\n');

    %else
    %load satSignal2.mat;
    %satSignal = savedSatSignal;
    %fprintf('-> Satellite signals  . . . . . . . . loaded;\n');
end %if(~skip) 

%% Generation of noise for each satellite =================================

%fprintf('Noise added\n')
if(~skip)
    fprintf('\t-Noise  . . . . . . . . . . . . . . added;\n');
end
%% Generation of multipath components =====================================
if(~skip)
    fprintf('\t-Multipath components . . . . . . . added;\n');
end
%% Gereneration of receiver signal ========================================
% call genRcvSignal to add all signals and save in approp. file
if(~skip)
    genRcvSignal(satSignal, satellites, settings, skip);
    fprintf('\t-Receiver signal  . . . . . . . . . created;\n');
    disp(['-> Signal Generator is over (elapsed time ', ...
                                        datestr(now - totalTime, 13), ')'])
    fprintf('\n\n');
else
    fprintf('-> Receiver signal  . . . . . . . . . loaded;\n');
end %if(skip)



