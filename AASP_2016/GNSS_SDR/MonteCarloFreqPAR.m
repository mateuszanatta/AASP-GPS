% MONTE CARLO SIMULATION FOR FREQUENCY ESTIMATION ERROR IN PARALLEL
% PROCESSSING MODE

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
settings.nrMSgen = 6;
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

SAT = [7 45 17653 0 -20 1 1];
    satellites = struct([]);
    satellites(1).PRN = SAT(1);
    satellites(1).DoA = SAT(2);
    satellites(1).CodPhase = SAT(3);
    satellites(1).FreqOffSet = SAT(4);
    satellites(1).SNR = SAT(5);
    satellites(1).Mpath = SAT(6);
    satellites(1).SIR = SAT(7);
    
    satCAtable = genSatCAtable(satellites);
    satNAVtable = genSatNAVtable(satellites);
    satPtable = genSatPtable();
    
    fprintf('\n')
    disp (['Monte Carlo Simulation started at ', datestr(totalTime)]);
    
    poolobj = gcp;
    addAttachedFiles(poolobj,{'satelliteFunctions','include'})
    
for bandFreq = 1:searchBinLength
    startTime = now;
    nowFreq = searchBin(bandFreq);
    satellites(1).FreqOffSet = nowFreq;
    fprintf('\n    Round %2d of %d frequencies\n',bandFreq,searchBinLength)
    disp (['    Generation started at ', datestr(startTime)]);
    skip = 0;
    
    %======================================================================
    %generateSignalMC
    samplesPerCode = round(settings.samplingFreq*settings.nrMSgen...
        /(settings.codeFreqBasis / settings.codeLength));
    satSignal = zeros(1,samplesPerCode);
    fprintf('Sat %2d ', satellites(1).PRN)
    satSignal(1,:) = genSatSignal(satellites(1), settings, samplesPerCode,...
            satCAtable(1,:), satNAVtable(1,:), satPtable);
    fprintf(' done;\n')
    noisedSignal = addNoise(satSignal,satellites,settings);
    data = genRcvSignalMC(noisedSignal, satellites, settings);
    %======================================================================
    
    IF = settings.IF;
    Phase = satellites.CodPhase;
    
    fprintf('   - Iteration 1 of %d \n',MCiterations)
    
    %postShortProcMC_PAR
    %======================================================================        
    samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
    longSignal = data;
    signal1 = longSignal(1 : samplesPerCode);
    signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);

    signal0DC = longSignal - mean(longSignal); 
    ts = 1 / settings.samplingFreq;
    phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;
    numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;
    caCodesTable = makeCaTable(settings);
    results     = zeros(numberOfFrqBins, samplesPerCode);
    frqBins     = zeros(1, numberOfFrqBins);
    for PRN = settings.acqSatelliteList
        caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));
        
        %==================================================================
        for frqBinIndex = 1:numberOfFrqBins
        %==================================================================
            frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               0.5e3 * (frqBinIndex - 1);
            sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
            cosCarr = cos(frqBins(frqBinIndex) * phasePoints);
            I1      = sinCarr .* signal1;
            Q1      = cosCarr .* signal1;
            I2      = sinCarr .* signal2;
            Q2      = cosCarr .* signal2;
            IQfreqDom1 = fft(I1 + 1j*Q1);
            IQfreqDom2 = fft(I2 + 1j*Q2);
            convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
            convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;
            acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
            acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;

            if (max(acqRes1) > max(acqRes2))
                results(frqBinIndex, :) = acqRes1;
            else
                results(frqBinIndex, :) = acqRes2;
            end    
        end % frqBinIndex = 1:numberOfFrqBins

        [~, frequencyBinIndex] = max(max(results, [], 2));
        [peakSize, codePhase] = max(max(results));
        samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
        excludeRangeIndex1 = codePhase - samplesPerCodeChip;
        excludeRangeIndex2 = codePhase + samplesPerCodeChip;
        if excludeRangeIndex1 < 2
            if (excludeRangeIndex1 == 1 || excludeRangeIndex1 == 0)
                 codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode - excludeRangeIndex1);
            else 
                 codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
            end
                         
        elseif excludeRangeIndex2 >= samplesPerCode
            codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
        else
            codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
        end

    
        secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));
        acqResults = struct;
        acqResults.peakMetric = peakSize/secondPeakSize;
    
        if (peakSize/secondPeakSize) > settings.acqThreshold
            caCode = generateCAcode(PRN);
            codeValueIndex = floor((ts * (1:5*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
            longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
            xCarrier = ...
              signal0DC(codePhase:(codePhase + 5*samplesPerCode-1)) ...
                .* longCaCode;
            fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
            fftxc = abs(fft(xCarrier, fftNumPts)); 
            uniqFftPts = ceil((fftNumPts + 1) / 2);
                [fftMax, fftMaxIndex] = max(fftxc(5 : uniqFftPts-5));
            fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
            acqResults.carrFreq  = fftFreqBins(fftMaxIndex);
            acqResults.codePhase = codePhase;
        else
            acqResults.carrFreq = 0;
            acqResults.codePhase = 0;
        end   % if (peakSize/secondPeakSize) > settings.acqThreshold  
    end    % for PRN = satelliteList

% =========================================================================
       
    errorFreq(bandFreq,1) = (acqResults.carrFreq - IF) - ...
                nowFreq;
    errorCode(bandFreq,1) = (acqResults.codePhase - ...
            Phase);
    clear('acqResults');
    
    parfor MC = 2:MCiterations
    
        fprintf('   - Iteration %2d of %d \n',MC,MCiterations)
        skip = 1;
        noisedSignal = addNoise(satSignal,satellites,settings);
        data = genRcvSignalMC(noisedSignal, satellites, settings);
        

        samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
        
        %postShortProcMC_PAR
        % =========================================================================        

        longSignal = data;
        signal1 = longSignal(1 : samplesPerCode);
        signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);
        signal0DC = longSignal - mean(longSignal); 
        ts = 1 / settings.samplingFreq;
        phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;
        numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;
        caCodesTable = makeCaTable(settings);
        results     = zeros(numberOfFrqBins, samplesPerCode);
        frqBins     = zeros(1, numberOfFrqBins);
        for PRN = settings.acqSatelliteList
            caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));

            for frqBinIndex = 1:numberOfFrqBins
                frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               0.5e3 * (frqBinIndex - 1);
                sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
                cosCarr = cos(frqBins(frqBinIndex) * phasePoints);
                I1      = sinCarr .* signal1;
                Q1      = cosCarr .* signal1;
                I2      = sinCarr .* signal2;
                Q2      = cosCarr .* signal2;
                IQfreqDom1 = fft(I1 + 1j*Q1);
                IQfreqDom2 = fft(I2 + 1j*Q2);
                convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
                convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;
                acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
                acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
                if (max(acqRes1) > max(acqRes2))
                    results(frqBinIndex, :) = acqRes1;
                else
                    results(frqBinIndex, :) = acqRes2;
                end
            end % frqBinIndex = 1:numberOfFrqBins

            [~, frequencyBinIndex] = max(max(results, [], 2));
            [peakSize, codePhase] = max(max(results));
            samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
            excludeRangeIndex1 = codePhase - samplesPerCodeChip;
            excludeRangeIndex2 = codePhase + samplesPerCodeChip;
        
            if excludeRangeIndex1 < 2
                if (excludeRangeIndex1 == 1 || excludeRangeIndex1 == 0)
                    codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode - excludeRangeIndex1);
                else 
                    codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
                end     
            elseif excludeRangeIndex2 >= samplesPerCode
                codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
            else
                codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
            end
            secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));
            acqResults = struct;
            acqResults.peakMetric = peakSize/secondPeakSize;
    
            if (peakSize/secondPeakSize) > settings.acqThreshold
                caCode = generateCAcode(PRN);
                codeValueIndex = floor((ts * (1:5*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
                longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
    
                xCarrier = ...
                    signal0DC(codePhase:(codePhase + 5*samplesPerCode-1)) ...
                    .* longCaCode;
                fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
                fftxc = abs(fft(xCarrier, fftNumPts)); 
                uniqFftPts = ceil((fftNumPts + 1) / 2);
                [fftMax, fftMaxIndex] = max(fftxc(5 : uniqFftPts-5));
                fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
                acqResults.carrFreq  = fftFreqBins(fftMaxIndex);
                acqResults.codePhase = codePhase;
            else
                acqResults.carrFreq = 0;
                acqResults.codePhase = 0;
            end   % if (peakSize/secondPeakSize) > settings.acqThreshold
        end    % for PRN = satelliteList

    % =========================================================================
        %postShortProcMC_PAR
        
        errorFreq(bandFreq,MC) = (acqResults.carrFreq - IF) - ...
                nowFreq;
        errorCode(bandFreq,MC) = (acqResults.codePhase - ...
            Phase);
    end %parfor MC = 2:MCiterations
    
    fprintf('\n')
    disp(['    -> Frequency simulation is over (elapsed time ', ...
                                        datestr(now - startTime, 13), ')'])
    disp(['    -> Monte Carlo Simulation elapsed time ', ...
                                        datestr(now - totalTime, 13), ')'])                                
    fprintf('\n----------------------------------------------------------')
end %bandFreq = 1:searchBinLength




save('errorFreq.mat','errorFreq')
save('errorCode.mat','errorCode')
save('AllSim16_12_06')
fprintf('\n\n')
disp(['-> Monte Carlo Simulation is over (elapsed time ', ...
                                        datestr(now - totalTime, 13), ')'])
disp('===================================================================')