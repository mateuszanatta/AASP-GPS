function [mDelayDiff] = delay_difference(pEstimatedDelays, pDelayTx, pSNR, pSatAmount)
%delay_difference Calculate the differences of the estimated delays minus de
%delays used for transmitting the signals. Then, squares the result.
%% Input
%   - pEstimatedDelays - The estimated delays of each signal for each SNR
%                      value
%   - pDelays          - The Delay values used to generate the transmitted
%                      signal
%   -pSNR              - The Signal Noise Ratio used to generate the
%                      signals
%   -pSatAmount        - The number of satellite transmitting
%% Output
%   - mDelayDiff - The outcome of the stimated delays minus the real delays

%%
    %Initiate matrices
    [m, n] = size(pEstimatedDelays);
    mDelayDiff = zeros(m, n, pSatAmount, length(pSNR));
    %Calculate the squared difference beteween the estimated delays and delays used
    %for generating the signals
    for snr = 1:length(pSNR)
        for sat = 1:pSatAmount
            mDelayDiff(:, :, sat, snr) = (bsxfun(@minus, real(pEstimatedDelays), pDelayTx')).^2;
        end
    end
end

