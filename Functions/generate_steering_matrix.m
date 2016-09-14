function [A] = generate_steering_matrix(pPhase, pNum_Sensor)
%GENERATE_STEERING_MATRIX Generate a Steering Matrix for a Uniform Linear
%Array (ULA)
%% Input
%   pPhase       - Signal Phase
%   pNum_Sensors - Number of sensors in the ULA
%% Output
%   A - The steering matrix

%%
d = size(pPhase,2);
%
A = zeros(pNum_Sensor, d);
for i = 1:d
    A(:, i) = exp(1j*pPhase(:,i)*(((-pNum_Sensor+1)/2:(pNum_Sensor-1)/2)'));
end

end

