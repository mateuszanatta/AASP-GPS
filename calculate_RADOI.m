function [ d ] = calculate_RADOI( X )
    % Determine signal dimensions
    [M,N] = size(X);
    
    % Covariance Matrix
    Rhat = 1/N*(X)*(X');
    
    % EVD
    [~, eigenvalues] = eig(Rhat, 'vector');
    eigenvalues = sort(eigenvalues, 'descend');

    % Determine cost mi
    for kk = 1:M-1;
        mi(kk) = 1/(M-kk) * sum( eigenvalues(kk+1:end) );
    end
    
    % Determine cost alpha
    alpha = max( (transpose(eigenvalues(1:end-1)) - mi)./mi );
    alpha = 1/alpha;
    
    % Determine cost epslon and the RADOI function
    for kk = 1:M-1;
        epslon(kk) = 1 - alpha*(eigenvalues(kk)-mi(kk))/mi(kk);
    end
    
    % Determine eigenvalues constant from 2nd till last value
    eig_sum = sum(eigenvalues(2:end));
    eig_sum = 1/eig_sum;
    
    % Determine epslon constant taking all elements (from 1 to M-1)
    eps_sum = sum(epslon(1:end));
    eps_sum = 1/eps_sum;
    
    for kk = 1:M-1;
        RADOI(kk) = eigenvalues(kk+1)*eig_sum - epslon(kk)*eps_sum;
    end
    
    % Model order
    [~, d] = min(RADOI);
end

