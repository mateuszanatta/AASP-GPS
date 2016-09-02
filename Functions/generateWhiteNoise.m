function [white_Noise] = generateWhiteNoise(pSignalSize, pSNR)

%   Generate White Noise

%% Input
%   pSignalSize     - The size of the MxN Matrix that contains the signal 
%   pSNR            - The signal to noise ratio

%% Output
%   white_Noise     - Returns a vector/matrix/tensor with same size as
%                   the signal
%%    
    sSigma = 10^(-pSNR/20)/sqrt(2);

    % Draw ZMCSCG (zero mean) WHITE noise with variance sigma^2
    white_Noise = (randn(pSignalSize)+ 1*1i*randn(pSignalSize))*sSigma;
end