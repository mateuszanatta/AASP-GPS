clear
clc

% Signal Properties
c = 4; % constellation size

% Model dimensions (d sources, M sensors, N samples)
d = 2; % sources (model order)
M = 3; % sensors (min d)
N = 2^10; % samples (M << N)

% Noise characterics
SNR = 20; % SNR (dB)
r = 10^(-SNR/20); % Ratio between noise and signal
n = (r * randn(M, N) + 1i* r * randn(M, N))/sqrt(2); % Complex gaussian noise

% Baseband signal
C = QAMconst(c); % 4 "QAM" constellation
S = C(randi(c,d,N)); % baseband signals

% Mixing matrix
[A,theta] = VMM(M, d); % Vandermonde Mixing Matrix
X = A*S + n; % Mixture

% ESPRIT
[theta_hat, a] = ESPRIT(X,d);

% Rebuild mixing matrix
W = pinv(a); % Demixing matrix
