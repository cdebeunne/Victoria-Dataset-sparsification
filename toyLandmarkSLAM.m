%% Create a toy example

% We create a toy pose graph to understand basic concept
%
%              > l3     > l4
%            /    ^   /    
%           /     |  /      
%    |--> x1 ---> x2
%
% The observation model is z = R' * (l - x)

toy_pg = poseGraph;
meas_inf_proprio = [1.5 0 0 4 0 400];
meas_inf_extero = [10 0 10];

% Relative pose
addRelativePose(toy_pg, [1, 0, 0], meas_inf_proprio, 1, 2);

% Landmark measurement
addPointLandmark(toy_pg, [0.5, 0.5], meas_inf_extero, 1, 3);
addPointLandmark(toy_pg, [-0.5, 0.5], meas_inf_extero, 2, 3);
addPointLandmark(toy_pg, [0, 0.5], meas_inf_extero, 2, 4);

% TO IMPROVE
% add a map to the type of node to build an information matrix with
% different node sizes
pose_nodes = [1 2];
n_pose = 2;
lmk_nodes = [3 4];
n_lmk = 2;

figure
show(toy_pg,'IDs','off');
title('toy graph');

%% Choose the node we want to remove and bookeping

node_to_remove = 2;
nodes_to_keep = [1 3 4];

%% Marginalization

% get the global markov blanket 
markov_blanket_remove = [];
markov_blanket_remove = [markov_blanket_remove,...
                       getMarkovBlanket(toy_pg, node_to_remove)];
markov_blanket_remove = unique(markov_blanket_remove);
n_mb = length(markov_blanket_remove);

% compute the information matrix on the elimination clique
I_t = computeMatInfJac(toy_pg, markov_blanket_remove);

% select the nodes to keep in the markov blanket
nodes_to_keep_t = nodes_to_keep(ismember(nodes_to_keep,markov_blanket_remove));
n_keep_t = length(nodes_to_keep_t);

% We select the rows and the columns to keep and to remove
rows_to_keep = [];
rows_to_remove = [];
for k = 1:length(nodes_to_keep_t)
    j = nodes_to_keep_t(k);
    if (ismember(j, pose_nodes))
        rows_to_keep = [rows_to_keep, (j-1)*3+1, (j-1)*3+2,...
            (j-1)*3+3];
    else
        s_q = n_pose * 3 + (j-n_pose - 1) * 2 + 1;
        rows_to_keep = [rows_to_keep, s_q, s_q+1];
    end
end
for k = 1:length(node_to_remove)
    j = node_to_remove(k);
    if (ismember(j, pose_nodes))
        rows_to_remove = [rows_to_remove, (j-1)*3+1, (j-1)*3+2,...
            (j-1)*3+3];
    else
        % s_q is the index in the information matrix of the landmark
        s_q = n_pose * 3 + (j-n_pose - 1) * 2 + 1;
        rows_to_remove = [rows_to_remove, s_q, s_q+1];
    end
end

% Compute the marginalized Information matrix with Schur Complement
I_aa = I_t(rows_to_keep, rows_to_keep);
I_bb = I_t(rows_to_remove, rows_to_remove);
I_ab = I_t(rows_to_keep, rows_to_remove);

I_marg = I_aa - I_ab * pinv(I_bb) * I_ab';


%% Select Manually a topology for the new subgraph

%            > l4  
%          /        
%         /         
%    |--x1 
%         \          
%          \        
%            > l3  
% This is where there is room for research 
% BEWARE the closed form solution cannot be used with loops!

% We use the rank revealing decomposition from Mazuran et al. 
[U, D] = eig(I_marg);
non_zero_columns = real(diag(D)) > 1e-5;
D = D(non_zero_columns, non_zero_columns);
U = U(:, non_zero_columns);
Sigma = inv(D);

% we define the topology
edges_in_tree = [1 4; 1 3];

% new nodes position 
% 1 -> 1
% 3 -> 2
% 4 -> 3
new_nodes = [1, 3, 4];

% We build Omega and J and retrieve measurements
J = zeros(size(I_marg));
Omega = zeros(size(I_marg));
n_pose = 1;

% compute info mat only for landmark factors
for k = 1:length(edges_in_tree)

    node_pair = edges_in_tree(k,:);
    i = node_pair(1);
    j = node_pair(2);
    
    % retrieve the robot pose
    meas_i = nodeEstimates(toy_pg,i);
    s = sin(meas_i(3));
    c = cos(meas_i(3));
    x = meas_i(1:2);
    R = [c -s;
        s c];

    % retrieve the landmark pose
    meas_j = nodeEstimates(toy_pg,j);
    l = meas_j(1:2);

    % Retrieve indices in the information matrix
    p = find(new_nodes==i);
    q = find(new_nodes==j);

    % s_q is the index in the information matrix of the landmark
    s_q = n_pose * 3 + (q-n_pose - 1) * 2 + 1;

    % Compute 2D jacobians
    % Formula from debeunne et al that is probably wrong
    J_k = zeros(2, length(I_marg));
    dz_dtheta = -[-s * x(1) + c * x(2); -c * x(1) - s * x(2)] +...
        [-s * l(1) + c * l(2); -c * l(1) - s * l(2)];
    dz_dx = -R';
    J_k(:, (p-1)*3+1:p*3) = [dz_dx dz_dtheta];
    J_k(:, s_q:s_q+1) = R';

    J(s_q:s_q+1, :) = J_k;

end

% Add unary prior
J(1:3, 1:3) = eye(3,3);

% Project on rank revealed subspace
J = J * U;

% Compute new info mat
inter = J * Sigma * J';

% Fill in stacked covariance
Omega(1:3, 1:3) = inv(inter(1:3, 1:3));
Omega(4:5, 4:5) = inv(inter(4:5, 4:5));
Omega(6:7, 6:7) = inv(inter(6:7, 6:7));

% Now we can compute the information matrix of the sparsified elimination
% clique

I_spars = J' * Omega * J;
I_spars = U * I_spars * U';

% plot it
figure;
imagesc(I_spars > 0);
title("Sparsified Information matrix on elimination clique");

%% Compute KLD between sparse and dense distribution

prod = I_spars * Sigma;
D_kl = 0.5 * ( trace(prod) - log(det(prod)) - length(I_spars))


