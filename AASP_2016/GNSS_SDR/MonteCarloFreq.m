% MONTE CARLO SIMULATION FOR FREQUENCY ESTIMATION ERROR

clear; close all; clc;
format ('compact');
format ('long', 'g');
addpath include

%% Monte Carlo settings ===================================================
MCiterations = 100;
searchBand = 10e3;
searchBandStep = 10;
searchBin = (-searchBand:searchBandStep:searchBand);
searchBinLength = length(searchBin);

%% Satellite definitions ==================================================
settings = initSettings();
settings.nrMSgen = 11;
settings.fileName           = ...
   'GPSsignalMCsim';
settings.path = '..\SatelliteSignals\MonteCarlo_records\';
addpath 'satelliteFunctions'

choose = 3;
if(choose == 2)
    fileNameStr = cat(2,settings.path,settings.fileName);
else
    fileNameStr = cat(2,settings.path,settings.fileName,'.bin');
end %(choose == 2)

if ~ischar(fileNameStr)
    error('File name must be a string');
end

errorFreq = zeros(searchBinLength,MCiterations);
errorCode = zeros(searchBinLength,MCiterations);

totalTime = now;
[fid, message] = fopen(fileNameStr, 'rb');
if fid > 0
    fclose(fid);
    delete(fileNameStr);
end

    fprintf('\n')
    disp (['Monte Carlo Simulation started at ', datestr(totalTime)]);
       
for bandFreq = 1:searchBinLength
    
    startTime = now;
    nowFreq = searchBin(bandFreq);
    SAT = [7 45 17653 nowFreq -20 1 1];
    fprintf('\n    Round %d of %d frequencies\n',bandFreq,searchBinLength)
    disp (['    Generation started at ', datestr(startTime)]);
    
    fprintf('   - Iteration 1 of %d \n',MCiterations)
    skip = 0;
    generateSignalMC
    postShortProcMC
    errorFreq(bandFreq,1) = (channel(1).acquiredFreq - settings.IF) - ...
        nowFreq;
    errorCode(bandFreq,1) = (channel(1).codePhase - ...
        satellites.CodPhase);
%    delete(fileNameStr);
    %{
    for MC = 2:MCiterations
    
        fprintf('   - Iteration %d of %d \n',MC,MCiterations)
        generateSignalMC
        postShortProcMC
        errorFreq(bandFreq,MC) = (channel(1).acquiredFreq - settings.IF) - ...
                searchBin(bandFreq);
        errorCode(bandFreq,MC) = (channel(1).codePhase - ...
            satellites.CodPhase);
%        delete(fileNameStr);
    end %for MC = 2:MCiterations
    %}
    
    IF = settings.IF;
    poolobj = gcp;
    addAttachedFiles(poolobj,{'satelliteFunctions','include'})
    parfor MC = 2:MCiterations
    
        %fprintf('   - Iteration %d of %d \n',MC,MCiterations)
        skip = 1;
        generateSignalMC
        postShortProcMC
        errorFreq(bandFreq,MC) = (channel(1).acquiredFreq - IF) - ...
                nowFreq;
        errorCode(bandFreq,MC) = (channel(1).codePhase - ...
            satellites.CodPhase);
%        delete(fileNameStr);
    end %for MC = 2:MCiterations
    
    disp(['    -> Frequency simulation is over (elapsed time ', ...
                                        datestr(now - startTime, 13), ')'])
    fprintf('\n----------------------------------------------------------')
end %bandFreq = 1:searchBinLength



%{    
for bandFreq = 1:searchBinLength
    startTime = now;
    nowFreq = searchBin(bandFreq);
    SAT = [7 45 17653 nowFreq -20 1 1];
    fprintf('\n    Round %d of %d frequencies\n',bandFreq,searchBinLength)
    disp (['    Generation started at ', datestr(startTime)]);
    
    fprintf('   - Iteration 1 of %d \n',MCiterations)
    skip = 0;
    generateSignalMC
    postShortProcMC
    errorFreq(bandFreq,1) = (channel(1).acquiredFreq - settings.IF) - ...
        nowFreq;
    errorCode(bandFreq,1) = (channel(1).codePhase - ...
        satellites.CodPhase);
%    delete(fileNameStr);
    skip = 1;
    %{
    for MC = 2:MCiterations
    
        fprintf('   - Iteration %d of %d \n',MC,MCiterations)
        generateSignalMC
        postShortProcMC
        errorFreq(bandFreq,MC) = (channel(1).acquiredFreq - settings.IF) - ...
                searchBin(bandFreq);
        errorCode(bandFreq,MC) = (channel(1).codePhase - ...
            satellites.CodPhase);
%        delete(fileNameStr);
    end %for MC = 2:MCiterations
    %}
    
    IF = settings.IF;
    for MC = 2:MCiterations
    
        %fprintf('   - Iteration %d of %d \n',MC,MCiterations)
        generateSignalMC
        postShortProcMC
        errorFreq(bandFreq,MC) = (channel(1).acquiredFreq - IF) - ...
                nowFreq;
        errorCode(bandFreq,MC) = (channel(1).codePhase - ...
            satellites.CodPhase);
%        delete(fileNameStr);
        disp(['    -> Frequency simulation is over (elapsed time ', ...
                                        datestr(now - startTime, 13), ')'])
    end %for MC = 2:MCiterations
    fprintf('\n----------------------------------------------------------')
end %bandFreq = 1:searchBinLength
%}
save('errorFreq.mat','errorFreq')
save('errorCode.mat','errorCode')

disp(['-> Monte Carlo Simulation is over (elapsed time ', ...
                                        datestr(now - totalTime, 13), ')'])
disp('===================================================================')