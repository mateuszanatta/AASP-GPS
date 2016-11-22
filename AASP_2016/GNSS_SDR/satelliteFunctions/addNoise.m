function signal = addNoise( signal, satellites , settings)
%GENRCVSIGNAL 
% INPUTS:   signal        - the satellites generated signal with respectly noise;
%           sats          - satellite struct with their informations
%           set           - settings struct with save file path;

%% Verify Multipath Components ============================================

%% Add Noise for EACH Satellite ===========================================
% Different noise for diferent satellites, but same noise for multipath
sqrSignal = signal.^2;
[sats,samples] = size(signal);

    for indexSat = 1:sats
        %{
        %measure signal power
        signalPow = sum(sqrSignal(indexSat,:)) / samples ;
    
        %calculate the desired noise power regarding actual signal power
        desirNoisePow = ( signalPow / db2mag( satellites(indexSat).SNR ) ); %in dB
    
        %create noise
        N = randn(1,samples);
    
        %measure noise power
        sqrN = N.^2;
        noisePow = sum( sqrN(:) ) / samples ;
    
        %apply desired noise power
        noiseAmp = sqrt( desirNoisePow / noisePow );
        N = noiseAmp .* N;
        %}
 %===========================================       
        %measure signal power
        signalPow = sum(sqrSignal(indexSat,:)) / samples ;
        
        %create noise
        N = randn(1,samples);
    
        %measure noise power
        sqrN = N.^2;
        noisePow = sum( sqrN(:) ) / samples ;
        %calculate the desired noise power regarding actual signal power
        desirSignalPow = ( noisePow * db2mag( satellites(indexSat).SNR ) ); %in dB
    
        %apply desired noise power
        signalAtten = sqrt( desirSignalPow / signalPow );
        signal(indexSat,:) = signal(indexSat,:) .* signalAtten;
%==============================================        
        %for m = 1:satellite(indexSat).Mpath
        signal(indexSat,:) = signal(indexSat,:) + N;
        %end %for m = 1:satellite(indexSat).Mpath
        
    end %for indexSat = 1:sats

end %function [] = addNoise( signal, satellites , settings)

