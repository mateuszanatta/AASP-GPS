function Us = est_sigsubsp_classic(x,d,dofba,dorvt)

% EST_SIGSUBSP_CLASSIC   Estimate signal subspace using matrix approach
%
%  US = EST_SIGSUBSP_CLASSIC(X,d[,F[,R]]) estimates the signal subspace
%  from the measurement matrix through an SVD. The rank of the signal
%  subspace must be provided by the parameter d.
%  The parameter F controls whether or not to use forward-backward
%  averaging (optional, defaults to 0). For F=1, R controls whether
%  or not to use the real-valued transformation. For F=0, R is ignored.
%  The resulting basis for the signal subspace is given in the matrix
%  US of size SIZE(X,1) by d.


if nargin < 4
    dorvt = 0;
end
if nargin < 3
    dofba = 0;
end
if dofba
    x = [x, fliplr(flipud(conj(x)))];
    if dorvt
        x = LPR_Matrix(size(x,1))'*x*LPR_Matrix(size(x,2));
    end
end
[U, S] = svd(x,0);
Us = U(:,1:d);