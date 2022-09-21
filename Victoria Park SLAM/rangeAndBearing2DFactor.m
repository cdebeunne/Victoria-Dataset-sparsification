function [l, I] = rangeAndBearing2DFactor(observation, covariance)
%RANGEANDBEARING2DFACTOR Summary of this function goes here
%   Detailed explanation goes here

% compute measurement
r = observation(1);
theta = observation(2);
l = [r * cos(theta); r * sin(theta)];

% compute information matrix
J = [cos(theta), -r * sin(theta);
     sin(theta), r * cos(theta)];
I = inv(J * covariance * J');


end

