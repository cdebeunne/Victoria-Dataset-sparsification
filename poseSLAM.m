%% Load g2o data

% Load posegraph
pg = load_g2o('input_MITb_g2o.g2o');
disp(pg);
show(pg,'IDs','off');
title('Original Pose Graph');

% Solve posegraph
updated_pg = optimizePoseGraph(pg, 'builtin-trust-region',...
    'VerboseOutput', 'on',...
    'InitialTrustRegionRadius',100);
figure
show(updated_pg,'IDs','off');
title('Updated Pose Graph');

%% Select nodes to keep and to remove 

node_indices = (1:updated_pg.NumNodes);
nodes_to_remove = [];
nodes_to_keep = [];

% node pairs that close loopsss
lc_node_pairs = updated_pg.edgeNodePairs(updated_pg.LoopClosureEdgeIDs);

p = 0.1; % probability of removing a node
for i=1:updated_pg.NumNodes

    % we do not remove nodes with lc
    if sum(lc_node_pairs == i) > 0
        nodes_to_keep = [nodes_to_keep i];
        continue;
    end

    if(rand < p)
        nodes_to_remove = [nodes_to_remove i];
    else
        nodes_to_keep = [nodes_to_keep i];
    end

end

n_remove = length(nodes_to_remove);
n_keep = length(nodes_to_keep);

%% Select nodes (punctual policy)

node_indices = (1:updated_pg.NumNodes);
nodes_to_remove = [2, 3];
nodes_to_keep = setdiff(node_indices, nodes_to_remove);

n_remove = length(nodes_to_remove);
n_keep = length(nodes_to_keep);

%% Marginalize over the elimination clique

% get the global markov blanket 
markov_blanket_remove = [];
for k = 1:n_remove
    markov_blanket_remove = [markov_blanket_remove,...
                           getMarkovBlanket(updated_pg, nodes_to_remove(k))];
end
markov_blanket_remove = unique(markov_blanket_remove);
n_mb = length(markov_blanket_remove);

% here the ordering is the same as markov_blanket_keep
I_t = computeMatInfJac(updated_pg, markov_blanket_remove);

% build the maps
map_index_t = containers.Map(1:n_mb, markov_blanket_remove);
map_node_t = containers.Map(markov_blanket_remove, 1:n_mb);

% select the nodes to keep in the markov blanket
nodes_to_keep_t = nodes_to_keep(ismember(nodes_to_keep,markov_blanket_remove));
n_keep_t = length(nodes_to_keep_t);

% reorder of the Information matrix in this way: 
% [  keep,keep    keep,remove ;
% [ remove,keep  remove,remove]

for k = 1:n_keep_t
    % I copy ?
    i = map_node_t(nodes_to_keep_t(k));

    % exchange rows
    rows_k = I_t((k-1)*3+1:k*3, :);
    rows_i = I_t((i-1)*3+1:i*3, :);
    I_t((k-1)*3+1:k*3, :) = rows_i;
    I_t((i-1)*3+1:i*3, :) = rows_k;

    % exchange columns
    cols_k = I_t(:,(k-1)*3+1:k*3);
    cols_i = I_t(:, (i-1)*3+1:i*3);
    I_t(:,(k-1)*3+1:k*3) = cols_i;
    I_t(:, (i-1)*3+1:i*3) = cols_k;

    % update maps
    map_index_t(k) = i;
    map_index_t(i) = k;
    map_node_t(i) = k;
    map_node_t(k) = i;
end

% Compute the marginalized Information matrix with Schur Complement
I_aa = I_t(1:3*n_keep_t, 1:3*n_keep_t);
I_bb = I_t(3*n_keep_t+1:3*n_mb, 3*n_keep_t+1:3*n_mb);
I_ab = I_t(1:3*n_keep_t,3*n_keep_t+1:3*n_mb);

I_marg = I_aa - I_ab * pinv(I_bb) * I_ab';

% plot it
imagesc(I_marg > 0);

%% Compute Chow Liu Tree

n = length(nodes_to_keep_t);

% an array with each row like [i, j, MI(x_i, x_j)]
map_pair_MI = [];

% compute every mutual information pair
% we work in the order [keep, remove]
for k=length(nodes_to_keep_t):-1:1
    for j=1:k-1
        map_pair_MI= [map_pair_MI; [j k computeMutualInfo(I_marg, j, k)]];
    end
end

% sort pairs wrt MI
map_pair_MI_sorted = sortrows(map_pair_MI, 3, 'descend');

% build tree using Chow & Liu algorithm (1968)
edges_in_tree = [];
nodes_in_tree = [];

for k=1:size(map_pair_MI_sorted, 1)
    pair = map_pair_MI_sorted(k, 1:2);

    % if the next node doesn't make a loop, we had it in the tree
    if (~(ismember(pair(1), nodes_in_tree) && ismember(pair(2), nodes_in_tree)))
        edges_in_tree = [edges_in_tree; pair];
        nodes_in_tree = [nodes_in_tree pair];
    end

end

%% Factor recovery in closed form

% First let's compute the covariance matrix
Sigma = inv(I_marg);

% We build Omega and J and retrieve measurements
J = zeros(size(I_marg));
Omega = zeros(size(I_marg));
z = zeros(3, length(edges_in_tree));

for k = 1:size(edges_in_tree, 1)

    node_pair = edges_in_tree(k,:);
    i = node_pair(1);
    j = node_pair(2);

    % Retrieve node estimates
    meas_i = nodeEstimates(updated_pg, map_index_t(i));
    ti = meas_i(1:2)';
    Ri = [cos(meas_i(3)) -sin(meas_i(3));
        sin(meas_i(3)) cos(meas_i(3))];

    meas_j = nodeEstimates(updated_pg, map_index_t(j));
    tj = meas_j(1:2)';
    Rj = [cos(meas_j(3)) -sin(meas_j(3));
        sin(meas_j(3)) cos(meas_j(3))];
    Rperp = [0 1; -1 0];

    delta_t = Ri' * (tj - ti);
    delta_R = (Ri' * Rj);
    delta_theta = acos(delta_R(1,1));
    z(:, k) = [delta_t', delta_theta]';

    % Compute 2D jacobians
    % Formula from 2D poseSLAM in GTSAM Dellaert.
    J_k = zeros(3, length(I_marg));
    J_k(:, (i-1)*3+1:i*3) = [Rj'*Ri Rperp*Rj'*(ti-tj);
                                    0 0 1];
    J_k(:, (j-1)*3+1:j*3) = eye(3,3);

    Omega_k = inv(J_k * Sigma * J_k');

    Omega((k-1)*3+1:k*3, (k-1)*3+1:k*3) = Omega_k;
    J((k-1)*3+1:k*3, :) = J_k;

end 

% We need to add an absolute pose factor so that J is invertible and we can
% use the closed form solution

J_k = zeros(3, length(I_marg));
J_k(1:3, 1:3) = eye(3,3);
k = 3;

Omega_k = inv(J_k * Sigma * J_k');
Omega((k-1)*3+1:k*3, (k-1)*3+1:k*3) = Omega_k;
J((k-1)*3+1:k*3, :) = J_k;

% Now we can compute the information matrix of the sparsified elimination
% clique

I_spars = J' * Omega * J;

% plot it
figure;
imagesc(I_spars > 1e-3);
title("Sparsified Information matrix on elimination clique");

