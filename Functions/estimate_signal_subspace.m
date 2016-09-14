function [Us] = estimate_signal_subspace(X, d)
%ESTIMATE_SIGNAL_SUBSPACE Estimates de signal subspace of signal matrix
% X       - Signal matrix
% d       - Model Order


[U, ~] = svd(X);

%Signal subspace
Us = U(:, 1:d);

end

