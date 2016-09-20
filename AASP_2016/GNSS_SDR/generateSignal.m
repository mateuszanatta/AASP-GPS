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
        'columns:  PRN       DoA      CodePhase      DopplerError\n',...
        '        (1 - 32)  (0º-90º) (1 - 1023chips) (-10k ~ +10k)Hz'])
    
satList = input('\nMatrix = \n');
% creating struct with satellites informations
satellites = [];
[d,~] = size(satList);
for ii = 1:d
    satellites(ii).PRN = satList(ii,1);
    satellites(ii).DoA = satList(ii,2);
    satellites(ii).CodPhase = satList(ii,3);
    satellites(ii).DoppErr = satList(ii,4);
end
showSatStatus(satellites);
fprintf(' Processing data: . . .\n');

%% Generation of satellite signals ========================================
% Preallocation of satSignal array
samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
satSignal = (samplesPerCode);
% call genSatSignal for each satellit
for ii = 1:d
    satSignal = genSatSignal(satellites(ii), settings, satSignal,...
         samplesPerCode);
end
%% Gereneration of receiver signal ========================================
% call genRcvSignal to add all signals and noise and save in approp. file


%% Generation 10 Navigation Tables as data for Simulink ===================
% Don't have preambles or valid navigations bits yet
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

%% Setup of the desired satellites =======================================================
