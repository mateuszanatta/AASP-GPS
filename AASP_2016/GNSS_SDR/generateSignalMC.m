%% generateSignalMC.m -- starts generations of GPS data
% Autor: Alexandre Serio Buscher ==========


%% Settings ===============================================================
% at settings struct

%% Ask inputs =============================================================
if(~skip)
satList = SAT;

% creating struct with satellites informations
[d,~] = size(satList);
satellites = struct([]);
for ii = 1:d
    satellites(ii).PRN = satList(ii,1);
    satellites(ii).DoA = satList(ii,2);
    satellites(ii).CodPhase = satList(ii,3);
    satellites(ii).FreqOffSet = satList(ii,4);
    satellites(ii).SNR = satList(ii,5);
    satellites(ii).Mpath = satList(ii,6);
    satellites(ii).SIR = satList(ii,7);
end %for ii = 1:d
satellites = orderStruct(satellites);
end
%% Generation of CA look-up table =========================================
if(~skip)
    satCAtable = genSatCAtable(satellites);
    %fprintf('\t-CA look-up table . . . . . . . . . created;\n');
end
%% Generation of Navigation lookup table ==================================
if(~skip)
    satNAVtable = genSatNAVtable(satellites);
    %fprintf('\t-Navigational look-up table . . . . created;\n');
end
%% Generation of P_yCode - simulated as Appendix B ========================
if(~skip)
    satPtable = genSatPtable(); %not generated
    %fprintf('\t-P(y) look-up table . . . . . . . . created;\n');
end
%% Generation of satellite signals ========================================
% Preallocation of satSignal array
%samplesPerCode = round(settings.samplingFreq*settings.nrMSgen*...
%    settings.nyquistGapgen/(settings.codeFreqBasis / settings.codeLength));

samplesPerCode = round(settings.samplingFreq*settings.nrMSgen...
    /(settings.codeFreqBasis / settings.codeLength));

if(~skip)
   
    satSignal = zeros(d,samplesPerCode);
 
    %fprintf('\t\t("." = 250us; "|" = 1 ms; "§" = 1 comp. C/A code)\n')
    % call genSatSignal to perform the signal composition
    for ii = 1:d
        fprintf('Sat %2d ', ii)
        satSignal(ii,:) = genSatSignal(satellites(ii), settings, samplesPerCode,...
            satCAtable(ii,:), satNAVtable(ii,:), satPtable);
        fprintf(' done;\n')
    end % for ii = 1:d
    %save(cat(2,settings.path,'\satSignal.mat'),'satSignal');
    %disp(['        Generation is over (elapsed time ', ...
    %                                    datestr(now - startTime, 13), ')'])
    %fprintf('\t-Satellite signals  . . . . . . . . created;\n');

end %if(~skip) 


%% Generation of noise for each satellite =================================

%fprintf('Noise added\n')
%if(~skip)
%    fprintf('Add noise?')
%    noise = input('   1 = Yes    0 = No\n');
%    if (noise)
        noisedSignal = addNoise(satSignal,satellites,settings);
%    end %if(noise)
%    fprintf('\t-Noise  . . . . . . . . . . . . . . added;\n');
%end
%% Gereneration of receiver signal ========================================
% call genRcvSignal to add all signals and save in approp. file
%if(~skip)
    data = genRcvSignalMC(noisedSignal, satellites, settings);
%    fprintf('\t-Receiver signal  . . . . . . . . . created;\n');
    %memory cleaning
%    clear d; clear ii; clear samplesPerCode; clear satCAtable;
%    clear satNAVtable; clear satPtable; %clear satSignal;
    
    
%%    clear startTime; clear totalTime;
    



