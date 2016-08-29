function X_ss = spatialsmooth_meastensor(X,L,sp)

% X    = [M1, M2, ..., MR, N]
% X_ss = [Msub1, Msub2, ..., MsubR, N*prod(L)]
% Msubr = Mr - L(r) + 1

% [10,1]
% [5]
% [6,5]

S = size(X);
R = length(S) - 1;

if length(L) ~= R
    S = [S, ones(1,length(L)-length(S)+1)]; % Append singleton dimensions for smoothing.
    R = length(S) - 1;
%    error('The length of the smoothing parameter vector should be equal to the number of dimensions in the measurement tensor minus one.');
end

if nargin < 3
    sp = 0;
end

M = S(1:R);
% N = S(R+1);

% Msub = M - L + 1;

%X_ss = zeros([Msub, N*L]);
X_ss = [];

l = ones(size(L));
if sp
    wbh = waitbar(0,'Smoothing...');
end
for nl = 1:prod(L)
    %l
    X_ss = cat(R+1, X_ss, full(subtensor(X, l, L, M)));
    if nl < prod(L)
        l(1) = l(1) + 1;
        for nnl = 1:length(l)
            if l(nnl) > L(nnl)
                l(nnl) = 1;
                l(nnl+1) = l(nnl+1)+1;
                %else
               % break
            end
        end
    end
    if sp
        waitbar(nl/prod(L),wbh);
    end
end
if sp
    close(wbh);
end
    


function X_s = subtensor(X, l, L, M)

R = length(M);

X_s = X;
for n = 1:length(l)
    Msub_n = M(n)-L(n)+1;
    sel_l = spalloc(Msub_n,M(n),Msub_n);
    sel_l(:,l(n):l(n)+Msub_n-1) = speye(Msub_n);
    if prod(size(sel_l))>1
        X_s = nmode_product(X_s,sel_l,n);
    end
end