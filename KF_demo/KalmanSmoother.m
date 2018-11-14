function [ smoothedxi,smoothedP ] = KalmanFilter( parameter_vector,y )
% This function runs the Kalman filter/smoother for the scalar AR(1) model plus
% noise with a diffuse prior (roughly)

% Extract lenght of the data
T = size(y,2);

% Extract the stuff we need from the input arguments
F = parameter_vector(1,1);
Q = parameter_vector(2,1);
R = parameter_vector(3,1);

% Run the Kalman filter
[xi,P,predictedxi,predictedP]=KalmanFilter(parameter_vector,y);

% Run the Kalman smoother 
for j=0:(T-1)
    t = T-j;
    % Set the last smoothed estimates equal to the filtered estimates
    if t>T-1
    smoothedxi(t) = xi(T);
    smoothedP(t)  = P(T);
    else
    smoothedxi(t) = xi(t) + P(t) * F' * 1/predictedP(t+1) * ( smoothedxi(t+1) - predictedxi(t+1) ) ;
    smoothedP(t)  = P(t)  - P(t) * F' * inv( predictedP(t+1) ) * ( predictedP(t+1) - smoothedP(t+1) ) * inv( predictedP(t+1) ) * F * P(t) ;
    end
end
% Close the function
end

