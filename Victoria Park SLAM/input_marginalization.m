function [z, Omega] = input_marginalization(meas_last, meas_curr, meas_new, I_meas_last, I_meas_curr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% compute the information matrix on the elimination clique

% Reshape measurements
meas_curr = reshape(meas_curr, 1, 3);
meas_last = reshape(meas_last, 1, 3);
meas_new = reshape(meas_new, 1, 3);


I_mb = zeros(9, 9);

% First factor
t_last = meas_last(1:2)';
R_last = [cos(meas_last(3)) -sin(meas_last(3));
    sin(meas_last(3)) cos(meas_last(3))];

t_curr = meas_curr(1:2)';
R_curr = [cos(meas_curr(3)) -sin(meas_curr(3));
    sin(meas_curr(3)) cos(meas_curr(3))];
Rperp = [0 1; -1 0];
J = zeros(3, 9);
J(:, 1:3) = -[R_curr'*R_last Rperp*R_curr'*(t_last-t_curr);
              0 0 1];
J(:, 4:6) = eye(3,3);
I_mb = I_mb + J' * I_meas_last * J; 

% Second factor
t_new = meas_new(1:2)';
R_new = [cos(meas_new(3)) -sin(meas_new(3));
    sin(meas_new(3)) cos(meas_new(3))];
J = zeros(3, 9);
J(:, 4:6) = -[R_new'*R_curr Rperp*R_new'*(t_curr-t_new);
              0 0 1];
J(:, 7:9) = eye(3,3);
I_mb = I_mb + J' * I_meas_curr * J; 

% We select the rows and the columns to keep and to remove
rows_to_keep = [1 2 3 7 8 9];
rows_to_remove = [4 5 6];

% Compute the marginalized Information matrix with Schur Complement
I_aa = I_mb(rows_to_keep, rows_to_keep);
I_bb = I_mb(rows_to_remove, rows_to_remove);
I_ab = I_mb(rows_to_keep, rows_to_remove);

I_marg = I_aa - I_ab * inv(I_bb) * I_ab';

%% Factor recovery in closed form

% We use the rank revealing decomposition from Mazuran et al. 
[U, D] = eig(I_marg);
non_zero_columns = real(diag(D)) > 1e-5;
D = D(non_zero_columns, non_zero_columns);
U = U(:, non_zero_columns);
Sigma = inv(D);

delta_t = R_last' * (t_new - t_last);
delta_R = (R_last' * R_new);
delta_theta = acos(delta_R(1,1));
z = [delta_t', delta_theta]';

% Compute 2D jacobians
% Formula from 2D poseSLAM in GTSAM Dellaert.
J = zeros(3, length(I_marg));
J(:, 1:3) = -[R_new'*R_last Rperp*R_new'*(t_last-t_new);
                                0 0 1];
J(:, 4:6) = eye(3,3);
J = J * U;

Omega = inv(J * Sigma * J');

end

