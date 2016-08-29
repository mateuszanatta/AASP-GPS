function T = GetVIT(M,K)

if (abs(K)==1)
    T=eye(M);
    if (sign(K)<0) 
        T = rot90(T);
    end;
    return;
end

r = (K+1) / (K-1);
ri = (1/abs(r)) * exp(1i*angle(r));
b = ((1-r)^(1-M)) * ((1-r)/(1-ri)).^[0:M-1];

a1 = ri*ones(M,1);
a2 = r*ones(M-1,1);
a1(1) = a2(1);

R=toeplitz(a1,a2);
T=zeros(M,M);

for k=1:M
    p = poly(R(k,:));
    T(k,:) = b(k) * p(end:-1:1);
end

