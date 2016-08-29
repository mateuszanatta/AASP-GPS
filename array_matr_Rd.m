function [A,D] = array_matr_Rd(mu,M,lpr)

% ARRAY_MATR_RD   Array matrices for R-D hypercuboidal array
%
% A = ARRAY_MATR_RD(MU,M) returns the array steering matrix
% for an R-D hypercuboidal array. Here R = LENGTH(M) and the 
% array has M(r) sensors in the r-th mode, r=1..R.
% MU contains the spatial frequencies and is of size R x d,
% where d is the number of sources. Each element of MU should
% be in [-pi, pi].
% A = ARRAY_MATR_RD(MU,M,LPR) allows additional control over
% the phase center of the array. The value of LPR controls
% which element is chosen as the phase center:
%    LPR = 0: the first element in each dimension [Default].
%    LPR = 1: the center element in each dimension. The
%    resulting array steering matrx is left-Pi-real.
% [A,D] = ARRAY_MATR_RD(MU,M[,LPR]) additionally returns the
% matrix of partial derivatives of the array steering vectors
% with respect to the spatial frequencies.

if nargin < 3
    lpr = 0;
end

R = size(mu,1);
d = size(mu,2);
if length(M) ~= R
    error('M should be a vector of size R, mu should be R x d.');
end
Mtotal = prod(M);

if any(mu>pi) | any(mu<-pi)
    warning('Spatial frequencies should be in [-pi, pi].');
end

A = zeros(Mtotal, d);
for i = 1:d
    Ac = cell(1,R);
    for r = 1:R
        if lpr
            Ac{r} = exp(j*mu(r,i)*(((-M(r)+1)/2:(M(r)-1)/2)'));
         else
            Ac{r} = exp(j*mu(r,i)*((0:M(r)-1)'));
        end
    end
    A(:,i) = Nkron(Ac);
end
        
D = zeros(Mtotal, R*d);
if nargout > 1
    for dr = 1:R
        for i = 1:d
            Dc = cell(1,R);
            for r = 1:R
                if r == dr
                    if lpr
                        Dc{r} = j*(((-M(r)+1)/2:(M(r)-1)/2)') .* exp(j*mu(r,i)*(((-M(r)+1)/2:(M(r)-1)/2)'));
                    else
                        Dc{r} = j*(0:M(r)-1)' .* exp(j*mu(r,i)*((0:M(r)-1)'));
                    end
                else
                    if lpr
                        Dc{r} = exp(j*mu(r,i)*(((-M(r)+1)/2:(M(r)-1)/2)'));
                    else
                        Dc{r} = exp(j*mu(r,i)*((0:M(r)-1)'));
                    end
                end
            end
            D(:,i+(dr-1)*d) = Nkron(Dc);
        end
    end
end