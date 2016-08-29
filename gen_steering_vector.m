function [ steeringVector ] = gen_steering_vector( M1, M2, phi, theta, xpos, ypos, B )
    % Generates Steering vector with different distance between sensors
    steeringVector = [];
    lambda = 3e8/B;
    for ii = 1:M1;
        for jj = 1:M2;
            % Include distance between antennas for steering vector calculation
            steeringVector = [steeringVector ...
                            1j*2*pi/lambda*...
                            (sqrt( xpos(ii)^2 + ypos(jj)^2 )*cos(phi)*sin(theta) + ...
                            sqrt( xpos(ii)^2 + ypos(jj)^2 )*sin(phi)*sin(theta))];
        end
    end
    steeringVector = transpose(exp(steeringVector));
end

