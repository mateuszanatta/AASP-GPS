function [mDelays] = delay_estimation(pDeviance, pLOS_Delays, pSNR, pSatAmount)
%DELAY_ESTIMATION Estimate the delay of a signal 
%   This function searches for the index of the minimal value encoutered in
%   the array pDeviance.
%% Inputs
%   - pDeviance   - The estimated deviance of the LOS of a Signal
%   - pSNR        - The Signal to Noise Ratio vector used.
%   - pLOS_Delays - A set of possible LOS delays for different SNR.
%   - pSatAmount  - Number of satellites trasmitting

%% Ouputs
%   - mDelays     - The Delay retrieved

%%
    %Initiate the matrix
    mDelays = zeros(pSatAmount, length(pSNR));
    
    %Find the delays
    for snr = 1:length(pSNR)
        for sat = 1:pSatAmount
            %Find index of the minimum real value in the matrix. Thus, one gets the
            %index of the lowest deviance between the received signal and the
            %correlator bank
            [~, index] = min(real(pDeviance(:, :, sat, snr)));
            
            %Once the retrieved index refers to the lowest deviance,
            % we get the value of the possible LOS delay.
            mDelays(sat, snr) = pLOS_Delays(index);
        end
    end


end

