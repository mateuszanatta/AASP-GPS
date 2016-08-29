function T = build_ten_Rd(factors)

% BUILD_TEN_RD   Build an R-D PARAFAC tensor from its factors
%
% Syntax:
%    T = BUILD_TEN_RD(Factors) where Factors is a cell array of length R
%    containing the factor matrices for each of the modes r=1,2,...,R. The
%    size of each factor matrix should be [M(r),d] where d is the number of
%    components (i.e., the rank of the tensor) and T will be of size
%    [M(1),M(2),...,M(R)].
%
% Author:
%    Florian Roemer, Communications Resarch Lab, TU Ilmenau
% Date:
%    Dec 2007

R = length(factors);
if R == 3
    % faster for R=3
    T = iunfolding(factors{3}*krp(factors{1},factors{2}).',3,[size(factors{1},1),size(factors{2},1),size(factors{3},1)]);
elseif R == 4
    % faster for R=4    
    T = iunfolding(factors{4}*krp_Rd(factors(1:3)).',4,[size(factors{1},1),size(factors{2},1),size(factors{3},1),size(factors{4},1)]);
else
    M = zeros(1,R);
    [M(1),d] = size(factors{1});

    for r = 2:R
        if size(factors{r},2) ~= d
            error('The number of columns must agree for all the factors');
        end
        M(r) = size(factors{r},1);
    end

    T = supereye(R,d);
    for r = 1:R
        T = nmode_product(T,factors{r},r);
    end
end