function T_prime = GetFocusedVIT( theta, d, M, K )
%Get Focused VIT
%   Focus the VIT in a different angle then mu = 0
%   [ T_prime ] = GetFocusedVIT( mu, M, K )
%   T_prime = T_k * P_mu
%   theta    DOA azimuth
%   d        inter element distance in wave lengths
%   M        number of sensor elements
%   K        VIT phase amplification, opt = 0.8

% generate P
P = diag(exp(1i*2*pi*d*sin((theta*pi/180))*[0:M-1]));
T_prime = GetVIT(M,K)*P;

end