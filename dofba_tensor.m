function Z = dofba_tensor(X,dim)

% DOFBA_TENSOR   Apply forward-backward averaging to a measurement tensor
%
% Syntax:
%   Z = DOFBA_TENSOR(X,dim)
% 
% Input:
%   X - measurement tensor of size M_1 x M_2 x ... x M_R x N
%   dim - Number of dimensions in X. Optional, defaults to NDIMS(X)-1. For
%   an R-D problem, dim must be equal to R+1.
%
% Output:
%   Z - forward-backward averaged measurement tensor of 
%       size M_1 x M_2 x ... x M_R x 2N
%
%
% Note:
%    If N = 1, then dim needs to be specified otherwise DOFBA_TENSOR will
%    produce wrong results. Since singleton dimensions are ignored the
%    number of dimensions will otherwise be wrong.

S = size(X);
if nargin > 1
    if length(S) < dim
        S = [S, ones(1,dim - length(S))];
    end
end    
R = length(S)-1;
M = S(1:R);
N = S(R+1);


thesubs_r = cell(1,R+1);
for r = 1:R
    thesubs_r{r} = M(r):-1:1;
end
thesubs_r{R+1} = N:-1:1;
Z = cat(R+1,X,subsref(conj(X),struct('type','()','subs',{thesubs_r})));
    
% alternative: X = invunfold(fba_2d(unfold(X,3).'),3,M(1)) <- slower!