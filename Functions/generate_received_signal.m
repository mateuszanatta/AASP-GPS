function [mReceivedSignal] = generate_received_signal(pSignal, pNumAntennas, pSignalPhase)
% Generate the signal received by one antenna or a Uniform Linear Array (ULA)
%
% This function will receive as parameter a set of signals without noise, 
% the number of antennas, and the signal phase.
%% Input
%   pSignal        - The transmited signal without noise
%   pNumAntennas   - The Number of Antennas to receive the signal. When
%                  more than one antenna, it assumes the array format is a ULA
%   pSignalPhase   - The phase used to generate 

%% Output
%   mReceivedSignal - The signal after received by the antenna without
%                   noise

%%
    % Generate received signal for an unique antenna 
    if(nargin < 2)
        %When using an unique antenna model the signal is transmitted
        %through an additive channel
        mReceivedSignal = sum(pSignal,2); 
    else
        if(pNumAntennas >= 2)

            %Generate the array antenna model
            vSteeringMatrix = generate_steering_matrix(pSignalPhase, pNumAntennas);

            %Generate the received signal by the ULA - Without Noise
            mReceivedSignal = vSteeringMatrix*transpose(pSignal);
        else
            disp('The number of antennas has to be greater than one');
        end
    end
    
end

