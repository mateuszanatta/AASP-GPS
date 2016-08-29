function [ Z ] = forwardBackwardAveraging( X )
    % Determine dimensions and snapshots
    [M,N] = size(X);
    
    % Determine PI matrix M
    Pi_M = fliplr(eye(M));

    % Determine PI matrix N
    Pi_N = fliplr(eye(N));

    % Forward Backward Averaging Result
    Z = [X Pi_M*conj(X)*Pi_N];

end

