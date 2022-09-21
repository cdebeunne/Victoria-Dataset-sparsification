function [dt, I] = carInputFactor(controller_input, L, timestep, covariance)
%CARINPUTFACTOR Summary of this function goes here
%   Detailed explanation goes here

% turn into relative pose
phi = controller_input(2);
v = controller_input(1);
dt = [v * timestep * cos(phi);
      v * timestep * sin(phi);
      tan(phi) * v * timestep / L];

% compute information matrix
J =  [timestep * cos(phi), - v * timestep * sin(phi);
      timestep * sin(phi), v * timestep * cos(phi);
      tan(phi) * timestep / L, v * timestep / (cos(phi)^2 * L)];
cov_mat = J * covariance * J';

 
% covariance regularization
if (rank(cov_mat) < 3 )
    cov_mat = cov_mat + eye(3);
end

% compute information matrix
I = inv(cov_mat);

end

