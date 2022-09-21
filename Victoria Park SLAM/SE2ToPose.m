function pose = SE2ToPose(M)
%SE2TOPOSE Summary of this function goes here
%   Detailed explanation goes here
theta = atan2(M(2,1), M(1,1));
t = M(1:2, 3)';

pose = [t theta]';
end

