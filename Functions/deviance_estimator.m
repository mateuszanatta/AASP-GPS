function [mDeviance, mPhase] = deviance_estimator(pSatAmount, pSignal, ...
                                                  pSignalDelays, pEstimated_Angles,...
                                                  pCorrelation_Bank, pMultiAntenna)
%DEVIANCE_ESTIMATOR Estimates the signal deviance based on Maximun
%Likelihood Estimator (MLE). Optionaly estimates the phase for an one antenna
%receiver.
%% Input
%   pSatAmount         - Number of satellites
%   pSignal            - The signal to be estimated
%   pSignalDelays      - A set of possible delayed signals
%   pEstimated_Angles  - A set of possible angles of arrival
%   pCorrelation_Bank  - Correlation Bank of a set of possible LOS
%   pMultiAntenna      - Set whether it is being used array of antennas or
%                        a single antenna. (1 - Array, 0 - Single)

%% Output
%   mDeviance   - The estimated deviance of the LOS of a Signal
%   mPhase      - The estimated phase of a singal. Only used when
%               estimating the phase of a one antenna receiver

%%
    %Retrieve the number of possible LOS in the Correlation Bank
    [numberLOS, pNumDelays, ~] = size(pCorrelation_Bank);
    
    %Preallocate Matrices
    mCorrelated_Signal = zeros(1, pNumDelays, pSatAmount);
    if(pMultiAntenna == 0)
        mEstimated_Deviance_LOS = zeros(numberLOS, length(pEstimated_Angles), pSatAmount);
    else
        mEstimated_Deviance_LOS = zeros(numberLOS, 1, pSatAmount);
    end
    
    %Calculates the Maximum Likelihood Estimator
    for sat = 1:pSatAmount
        if(pMultiAntenna == 0)
            for vAngle = 1:length(pEstimated_Angles)
               %Correlate the signal with the Delayed signal 
               %Multiplies pSignal by pOffSetCorr in order to eliminate the
               %signal phase
                 mCorrelated_Signal(:, :, sat) = (pSignal * ...
                     exp(-1j*pEstimated_Angles(vAngle)))' * pSignalDelays(:, :, sat);

               %% Applies MLE in order to estimate signal delay
               %Takes the difference between the Correlation Bank (pCorrelationBank) and the
               %matrix of correlated signal (mCorrelated_Signal)
               corr_Banks_Diff = bsxfun(@minus, pCorrelation_Bank(:, :, sat), mCorrelated_Signal(:, :, sat));
               %corr_Banks_Diff = pCorrelation_Bank(:, :, sat) - replicated_Corr_Signal;

               %Estimates de deviance between the received signal and the set
               %of possible LOS in the correlation bank (pCorrelationBank)
%                mEstimated_Deviance_LOS(:, vAngle, sat) = sqrt(sum((corr_Banks_Diff.^2), 2)); 
               mEstimated_Deviance_LOS(:, vAngle, sat) = sum(sqrt(corr_Banks_Diff.^2),2);           

            end %END For
        else
            %Correlate the signal with the Delayed signal 
               %Multiplies pSignal by pOffSetCorr in order to eliminate the
               %signal phase
           mCorrelated_Signal(:, :, sat) = (pSignal(:,sat) * ... 
                     exp(-1j*pEstimated_Angles(sat)))' * pSignalDelays(:, :, sat);
             
           %% Applies MLE in order to estimate signal delay

           %Takes the difference between the Correlation Bank (pCorrelationBank) and the
           %matrix of correlated signal (mCorrelated_Signal)
           corr_Banks_Diff = bsxfun(@minus, pCorrelation_Bank(:, :, sat), mCorrelated_Signal(:,:,sat));
%            corr_Banks_Diff = pCorrelation_Bank(:, :, sat) - replicated_Corr_Signal;

           %Estimates de deviance between the received signal and the set
           %of possible LOS in the correlation bank (pCorrelationBank)
           mEstimated_Deviance_LOS(:, :, sat) = sqrt(sum((corr_Banks_Diff.^2), 2)); 
%            mEstimated_Deviance_LOS(:, :, sat) = sum(sqrt(corr_Banks_Diff.^2),2);
        end%END IF
        
    end %END For
    
    %% 
    %Finds the lowest value in the estimated deviance of LOS
    %The lowest value indicates the shortest distance between the received
    %signal and the correlation bank
    for sat = 1:pSatAmount
        %find the lowest value
        minDeviance = sum(sqrt(imag(mEstimated_Deviance_LOS(:, :, sat)).^2));
        %get the index of the lowest value
        [~, minDeviance_Index] = min(minDeviance);
        
        %store the estimated deviance to each satellite
        mDeviance(:, :, sat) = mEstimated_Deviance_LOS(:, minDeviance_Index, sat);
        
        %%
        % Phase estimation when using one antenna
        if nargout > 1
            % Estimating transmitted signal phases. This algorithm does not
            % consist in a classic DOA estimation algorithm, but only a Maximum
            % Likelihood Estimator to find the phase as well.
            mPhase(:,sat) = pEstimated_Angles(minDeviance_Index);
        end
    end %END For
end

