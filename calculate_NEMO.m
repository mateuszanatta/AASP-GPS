function [ d ] = calculate_NEMO( X )
    % Determine signal dimensions
    [M,N] = size(X);
    
    % Covariance Matrix
    Rhat = 1/N*(X)*(X');
    
    % EVD
    [~, eigenvalues] = eig(Rhat, 'vector');
    eigenvalues = sort(eigenvalues, 'descend');
    
    % Determine phi constant
    if ( isreal(X) )
        phi = 1;
    else
        phi = 2;
    end
    
    % Look for minimum argument in probability function
    for kk = 1:M-1
        t(kk) = ((M-kk) * (sum(eigenvalues(kk+1:end).^2) / sum(eigenvalues(kk+1:end))^2) - (1+M/N))*N - M/N*(2/phi -1);
        NEMO(kk) = phi/4*(N/M)*t(kk)^2 + 2*(kk+1);
    end
    
    % Model order
    [~, d] = min(NEMO);
end

