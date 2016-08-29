function [Factors,inno,recerr] = plainvanilla_tals_Rd(X,d,Factors0,MAXIT)

% PLAINVANILLA_TALS_RD   "Plain vanilla" implementation of multi-linear
% alternating least squares to solve for PARAFAC models.
%
% Syntax:
%    Factors = PLAINVANILLA_TALS_RD (X,d[,InitFactors[,MAXIT]])
% 
% Input:
%    X           - R-D data tensor
%    d           - number of factors
%    InitFactors - Initialization for the ALS. Optional, if not supplied, a
%    random initialization is used
%    MAXIT       - Max. number of iterations allowed. Optional, defaults to 10000.
% 
% Output:
%    Factors     - length-R cell array containing the R estimated PARAFAC
%       matrices. Factors{r} will be of size [M(r),d], where M = SIZE(X).
%
% Author:
%    Florian Roemer, Communications Resarch Lab, TU Ilmenau
% Date:
%    Dec 2007

M = size(X);
R = length(M);
if (nargin > 2) && (~isempty(Factors0))
    Factors = Factors0;
else
    Factors = cell(1,R);
    if isreal(X)
        for r = 1:R
            Factors{r} = randn(M(r),d);
        end
    else
        for r = 1:R
            Factors{r} = randn(M(r),d) + j*randn(M(r),d);
        end
    end
end

lastFactors = Factors;
if nargin < 4
    MAXIT = 10000;
end
delta = 1e-8;

Xu = unfoldings(X);
inno = zeros(1,MAXIT);
if nargout > 2
    recerr = zeros(1,MAXIT);
    Xn = ho_norm(X);
end
for nIt = 1:MAXIT
    for r = R:-1:1
        Factors{r} = Xu{r} * pinv(krp_Rd(Factors([r+1:R,1:r-1]))).';
    end
  
    for r = 1:R
        if any(isinf(Factors{r}(:)))
            inno(nIt) = -1;
            Factors = lastFactors;
        else
            inno(nIt) = inno(nIt) + norm(Factors{r}(:) - lastFactors{r}(:)) / norm(Factors{r}(:));
        end
    end
    if inno(nIt) < delta
        inno = inno(1:nIt);
        if nargout > 2
            recerr = recerr(1:nIt);
        end
        break;
    end
    lastFactors = Factors;
    if nargout > 2
        Xrec = build_ten_Rd(Factors);
        recerr(nIt) = ho_norm(X - Xrec) / Xn;
    end
end

