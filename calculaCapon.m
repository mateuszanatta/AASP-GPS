function [ Pw, theta, PwMusic ] = calculaCapon( X, Rxx, M, delta, lambda )
% Inversa da matriz de covariancia
Rxx_inv = inv(Rxx);

% Direcao de chegada
theta = -pi/2:1/(200*pi):pi/2;

for ii = 1:length(theta);
    % Steering vector em formato vandermonde
    a_sv(1:M,1) = exp(j*(0:M-1)*2*pi*delta*sin(theta(ii))/lambda);

    % Capon
    Pw(ii) = 1/(a_sv'*(Rxx_inv)*a_sv);
end


% SVD de X
[Usvd, ~, ~] = svd(X);
d = 3;
Un = Usvd(:,(d+1):end);

for ii = 1:length(theta);
    % Steering vector em formato vandermonde
    a_sv(1:M,1) = exp(j*(0:M-1)*2*pi*delta*sin(theta(ii))/lambda);

    % Capon
    PwMusic(ii) = 1/(a_sv'*(Un)*(Un')*a_sv);
end

end

