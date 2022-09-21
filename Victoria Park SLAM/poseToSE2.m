function M = poseToSE2(pose)
%POSETOSE2 Summary of this function goes here
%   Detailed explanation goes here

theta = pose(3);
R = [cos(theta), -sin(theta);
    sin(theta), cos(theta)];
t = [pose(1) pose(2)]';
M = [R t; 0 0 1];

end

