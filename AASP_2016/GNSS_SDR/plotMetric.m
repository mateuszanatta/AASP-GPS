function plotMetric(PRN,longSignal,settings)
%Functions plots bar plot of acquisition results (acquisition metrics). No
%bars are shown for the satellites not included in the acquisition list (in
%structure SETTINGS). 
%
%plotAcquisition(acqResults)
%
%   Inputs:
%       acqResults    - Acquisition results from function acquisition.
%       settings      - settings struct containing threshold values

%--------------------------------------------------------------------------

%% Recover metric information for desired satellite =======================
% As acquisition.m
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

% Create two 1msec vectors of data to correlate with and one with zero DC
signal1 = longSignal(1 : samplesPerCode);
signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);

signal0DC = longSignal - mean(longSignal); 

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = makeCaTable(settings);

%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, samplesPerCode);

% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);

caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));

    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins

        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
        frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               0.5e3 * (frqBinIndex - 1);

        %--- Generate local sine and cosine -------------------------------
        sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
        cosCarr = cos(frqBins(frqBinIndex) * phasePoints);

        %--- "Remove carrier" from the signal -----------------------------
        I1      = sinCarr .* signal1;
        Q1      = cosCarr .* signal1;
        I2      = sinCarr .* signal2;
        Q2      = cosCarr .* signal2;

        %--- Convert the baseband signal to frequency domain --------------
        IQfreqDom1 = fft(I1 + 1j*Q1);
        IQfreqDom2 = fft(I2 + 1j*Q2);

        %--- Multiplication in the frequency domain (correlation in time
        %domain)
        convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
        convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;

        %--- Perform inverse DFT and store correlation results ------------
        acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
        acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
        
        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
        if (max(acqRes1) > max(acqRes2))
            results(frqBinIndex, :) = acqRes1;
        else
            results(frqBinIndex, :) = acqRes2;
        end
        
    end % frqBinIndex = 1:numberOfFrqBins

%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    [~, frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize, codePhase] = max(max(results));

    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
       %orginal book's code BUG; add this
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

    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

    %--- Store result -----------------------------------------------------
    metric = peakSize/secondPeakSize;
    fprintf ('\n Satellite %2d metric result: %3.2f',PRN, metric)
    fprintf('\n')
    
    %% Fine resolution frequency search ===================================
caCode = generateCAcode(PRN);
        
        codeValueIndex = floor((ts * (1:10*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
                            
        longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
    
        %--- Remove C/A code modulation from the original signal ----------
        % (Using detected C/A code phase)
        xCarrier = ...
            signal0DC(codePhase:(codePhase + 10*samplesPerCode-1)) ...
            .* longCaCode;
        
        %--- Find the next highest power of two and increase by 8x --------
        fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
        
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency 
        fftxc = abs(fft(xCarrier, fftNumPts)); 
        
        uniqFftPts = ceil((fftNumPts + 1) / 2);
        [fftMax, fftMaxIndex] = max(fftxc(5 : uniqFftPts-5));
        
        fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
               
%% Plot all results =======================================================
figure(113);

h3D = meshz(results);
axes = gca;
title ('Correlation results')
axes.XLabel.String = 'Sample number';
axes.YLabel.String = {'Frequency', '[MHz]   '};
N = length(frqBins)-1;
L = length(axes.YTick);
P = N/(L-1);
yTickLabel = zeros(1,L);
yTickLabel(1) = frqBins(1);
for ii=2:L
    yTickLabel(ii) = frqBins( (ii-1)*P );
end
yTickLabel = yTickLabel*1e-6;
axes.YTickLabel = yTickLabel;
axes.ZLabel.String = 'Correlation magnitude';

valueFreq = num2str(frqBins(frequencyBinIndex)*1e-6);
text1 = cat(2,'Acq  Freq  : ',valueFreq,' MHz');
valuePhase = num2str(codePhase);
text2 = cat(2,'Acq Phase: ',valuePhase,' samp.');
annotation('textbox', [0.13,0.8,0.3,0.1],'String', {text1,text2})
%                [hrz_start vrt_start hrz_width vrt_width]

