function [ steeringVector ] = gen_steering_vector_simple( M1, M2, phi, B )
    % Generates Steering vector with different distance between sensors
    steeringVector = [];
    
    % Signal Wavelength
    lambda = 3e8/B;
    
    % Distance between antennas
    delta = lambda/2;
    
    for ii = 1:M1*M2;
        % Include distance between antennas for steering vector calculation
        steeringVector = [steeringVector ...
                          exp(1j*2*pi/lambda*delta*cos(phi)*(ii-1))];
    end
    steeringVector = transpose(steeringVector);
end