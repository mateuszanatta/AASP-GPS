function U = Nkron(UC)

% NKRON   N-fold Kronecker product
%
% Syntax:
%    U = NKRON(UC)
% 
% Input:
%    UC is expected as a cell array of length R, each element being a
%    matrix of arbitrary size.
%
% Output:
%    U will be a matrix computed by UC{1} \kron UC{2} \kron ... \kron UC{R}.

U = UC{1};
for n=2:length(UC)
    U = kron(U,UC{n});
end