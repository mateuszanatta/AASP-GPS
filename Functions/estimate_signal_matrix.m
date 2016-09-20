function [mEstimatedSignal] = estimate_signal_matrix(pSteeringMatrix, pMixingMatrix)
%ESTIMATE_SIGNAL_MATRIX Estimates the signal matrix
%   This function multiplicates the mixing matrix by the pseudo-inverse
%   of the Steering Matrix (A Matrix) 
%% Input
%   - pSteeringMatrix - Receives the Steering Matrix (A matrix) of the
%                     antenna array. This parameter can be an estimated
%                     steering matrix or the steering matrix used to create
%                     the antenna array.
%
%   - pSignal         - The received signal by the antenna
%% Output
%   - mEstimatedSignal - The received estimate received signal. This one is
%                      obtained by multiplying the pseudo-inverse of the 
%                      Steering Matrix by the received signal

    if(nargin >= 2)
        mEstimatedSignal = pinv(pSteeringMatrix)*pSignal;
    else
end

