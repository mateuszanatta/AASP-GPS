function [ d ] = calculate_AIC( X )
    % Received Signal dimensions and Snapshots
    [M,N] = size(X);
    
    % Estimated Covariance Matrix
    Rhat = 1/N*(X)*(X');
    
    % Eigenvalue decomposition
    [~, eigenvalues] = eig(Rhat, 'vector');
    eigenvalues = sort(eigenvalues,'descend');
    
    % Perform AIC estimator
    for ii = 1:length(eigenvalues);
        % numerator for likelihood function
        eigenvalues_num = eigenvalues(ii+1:end);
        prod_num = eigenvalues_num.^(1/(M-ii));
        prod_num = prod(prod_num);

        % denominator for likelihood function
        eigenvalues_den = eigenvalues(ii+1:end);
        sum_den = sum(eigenvalues_den);
        sum_den = (1/(M-ii))*sum_den;

        % likelihood function
        argAIC = (prod_num/sum_den)^(N*(M-ii));
        AIC(ii) = -2*log(argAIC) + 2*ii*(2*M-ii);
    end
    
    % Model Order
    [~, d] = min(AIC);
end

