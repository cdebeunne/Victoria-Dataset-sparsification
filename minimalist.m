%% Create a toy example

% We create a toy pose graph to understand basic concept
%
%                      > x3 
%                    /    
%                   /       
%        x1 ---> x2 

toy_pg = poseGraph;
measInf = [1.5 0 0 4 0 400];
I_meas = vec_to_mat(measInf);

addRelativePose(toy_pg, [1, 0, pi/4], measInf, 1, 2);
addRelativePose(toy_pg, [1, 0, -pi/4], measInf, 2, 3);

updated_toy_pg = optimizePoseGraph(toy_pg, 'builtin-trust-region',...
    'VerboseOutput', 'on',...
    'InitialTrustRegionRadius',100);

figure
show(updated_toy_pg,'IDs','off');
title('toy graph');

%% Choose the node we want to remove and bookeping

node_to_remove = 2;
nodes_to_keep = [1 3];

%% Marginalization on the elimination clique

% get the global markov blanket 
markov_blanket_remove = [];
markov_blanket_remove = [markov_blanket_remove,...
                       getMarkovBlanket(updated_toy_pg, node_to_remove)];
markov_blanket_remove = unique(markov_blanket_remove);
n_mb = length(markov_blanket_remove);

% compute the information matrix on the elimination clique
I_mb = zeros(9, 9);

% First factor
meas_last = nodeEstimates(updated_toy_pg, 1);
t_last = meas_last(1:2)';
R_last = [cos(meas_last(3)) -sin(meas_last(3));
    sin(meas_last(3)) cos(meas_last(3))];

meas_curr = nodeEstimates(updated_toy_pg, 2);
t_curr = meas_curr(1:2)';
R_curr = [cos(meas_curr(3)) -sin(meas_curr(3));
    sin(meas_curr(3)) cos(meas_curr(3))];
Rperp = [0 1; -1 0];
J = zeros(3, 9);
J(:, 1:3) = -[R_curr'*R_last Rperp*R_curr'*(t_last-t_curr);
              0 0 1];
J(:, 4:6) = eye(3,3);
I_mb = I_mb + J' * I_meas * J; 

% Second factor
meas_new = nodeEstimates(updated_toy_pg, 3);
t_new = meas_new(1:2)';
R_new = [cos(meas_new(3)) -sin(meas_new(3));
    sin(meas_new(3)) cos(meas_new(3))];
J = zeros(3, 9);
J(:, 4:6) = -[R_new'*R_curr Rperp*R_new'*(t_curr-t_new);
              0 0 1];
J(:, 7:9) = eye(3,3);
I_mb = I_mb + J' * I_meas * J; 

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

% Now we can compute the information matrix of the sparsified elimination
% clique


I_spars = J' * Omega * J;

% plot it
figure;
clims = [-6, 6];

subplot 121;
imagesc(log(abs(I_marg)), clims);
title("Marginalized, Dense Information matrix on elimination clique")

subplot 122;
imagesc(log(abs(I_spars)), clims);
title("Sparsified Information matrix on elimination clique");

%% test function

[z_test, Omega_test] = input_marginalization(meas_last, meas_curr, meas_new, I_meas);
