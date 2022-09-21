%% Create a toy example

% We create a toy pose graph to understand basic concept
%
%                      > x4 ---> x5
%                    /    ^
%                   /     |  
%        x1 ---> x2 --->  x3

toy_pg = poseGraph;
measInf = [1.5 0 0 4 0 400];

addRelativePose(toy_pg, [1, 0, pi/4], measInf, 1, 2);
addRelativePose(toy_pg, [sqrt(0.5), -sqrt(0.5), -pi/4], measInf, 2, 3);
addRelativePose(toy_pg, [sqrt(2), 0, -pi/4], measInf, 2, 4);
addRelativePose(toy_pg, [0, 1, 0], measInf, 3, 4);
addRelativePose(toy_pg, [1, 0, 0], measInf, 4, 5);

updated_toy_pg = optimizePoseGraph(toy_pg, 'builtin-trust-region',...
    'VerboseOutput', 'on',...
    'InitialTrustRegionRadius',100);

figure
show(updated_toy_pg,'IDs','off');
title('toy graph');

%% Choose the node we want to remove and bookeping

node_to_remove = 2;
nodes_to_keep = [1 3 4 5];

%% Marginalization on the elimination clique

% get the global markov blanket 
markov_blanket_remove = [];
markov_blanket_remove = [markov_blanket_remove,...
                       getMarkovBlanket(updated_toy_pg, node_to_remove)];
markov_blanket_remove = unique(markov_blanket_remove);
n_mb = length(markov_blanket_remove);

% compute the information matrix on the elimination clique
I_mb = computeMatInfJac(updated_toy_pg, markov_blanket_remove);

% select the nodes to keep in the markov blanket
nodes_to_keep_mb = nodes_to_keep(ismember(nodes_to_keep,markov_blanket_remove));
n_keep_t = length(nodes_to_keep_mb);

% We select the rows and the columns to keep and to remove
rows_to_keep = [];
rows_to_remove = [];
for k = 1:length(nodes_to_keep_mb)
    j = nodes_to_keep_mb(k);
    rows_to_keep = [rows_to_keep, (j-1)*3+1, (j-1)*3+2,...
        (j-1)*3+3];
end
for k = 1:length(node_to_remove)
    j = node_to_remove(k);
    rows_to_remove = [rows_to_remove, (j-1)*3+1, (j-1)*3+2,...
        (j-1)*3+3];
end

% Compute the marginalized Information matrix with Schur Complement
I_aa = I_mb(rows_to_keep, rows_to_keep);
I_bb = I_mb(rows_to_remove, rows_to_remove);
I_ab = I_mb(rows_to_keep, rows_to_remove);

I_marg = I_aa - I_ab * inv(I_bb) * I_ab';

%% Compute Chow Liu Tree

n = length(nodes_to_keep_mb);

% an array with each row like [i, j, MI(x_i, x_j)]
map_pair_MI = [];

% compute every mutual information pair
for k=length(nodes_to_keep_mb):-1:1
    for j=1:k-1
        map_pair_MI= [map_pair_MI; [j k computeMutualInfo(I_marg, j, k)]];
    end
end

% sort pairs wrt MI
map_pair_MI_sorted = sortrows(map_pair_MI, 3, 'descend');

% build tree using Chow & Liu algorithm (1968)
edges_in_tree = [];
nodes_in_tree = [];

for k=1:length(map_pair_MI_sorted)
    pair = map_pair_MI_sorted(k, 1:2);

    % if the next node doesn't make a loop, we had it in the tree
    if (~(ismember(pair(1), nodes_in_tree) && ismember(pair(2), nodes_in_tree)))
        edges_in_tree = [edges_in_tree; pair];
        nodes_in_tree = [nodes_in_tree pair];
    end

end

%% Factor recovery in closed form

% We use the rank revealing decomposition from Mazuran et al. 
[U, D] = eig(I_marg);
non_zero_columns = real(diag(D)) > 1e-5;
D = D(non_zero_columns, non_zero_columns);
U = U(:, non_zero_columns);
Sigma = inv(D);

% We build Omega and J and retrieve measurements
J = zeros(3*size(edges_in_tree, 2), length(I_marg));
Omega = zeros(3*size(edges_in_tree, 2), 3*size(edges_in_tree, 2));
z = zeros(3, length(edges_in_tree));

for k = 1:length(edges_in_tree)

    node_pair = edges_in_tree(k,:);
    i = node_pair(1);
    j = node_pair(2);

    % Retrieve node estimates
    meas_i = nodeEstimates(updated_toy_pg, nodes_to_keep_mb(i));
    ti = meas_i(1:2)';
    Ri = [cos(meas_i(3)) -sin(meas_i(3));
        sin(meas_i(3)) cos(meas_i(3))];

    meas_j = nodeEstimates(updated_toy_pg, nodes_to_keep_mb(j));
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
    J_k(:, (i-1)*3+1:i*3) = -[Rj'*Ri Rperp*Rj'*(ti-tj);
                                    0 0 1];
    J_k(:, (j-1)*3+1:j*3) = eye(3,3);


    J((k-1)*3+1:k*3, :) = J_k;

end 

J = J * U;

for k = 1:length(edges_in_tree)
    inter = J * Sigma * J';
    Omega_k = inv(inter((k-1)*3+1:k*3, (k-1)*3+1:k*3));
    Omega((k-1)*3+1:k*3, (k-1)*3+1:k*3) = Omega_k;
end


% Now we can compute the information matrix of the sparsified elimination
% clique

I_spars = J' * Omega * J;
I_spars = U * I_spars * U';

% plot it
figure;
clims = [-6, 6];

subplot 121;
imagesc(log(abs(I_marg)), clims);
title("Marginalized, Dense Information matrix on elimination clique")

subplot 122;
imagesc(log(abs(I_spars)), clims);
title("Sparsified Information matrix on elimination clique");


%% Compute KLD between sparse and dense distribution

prod = I_spars * Sigma;
D_kl = 0.5 * ( trace(prod) - log(det(prod)) - length(I_spars))

%% Build a new pose graph 

sparse_pg = poseGraph;

% we add the edges from the sparsified markov blanket
for k = length(edges_in_tree):-1:1
    addRelativePose(sparse_pg, z(:,k), mat_to_vec(Omega((k-1)*3+1:k*3, (k-1)*3+1:k*3)),...
        edges_in_tree(k,1), edges_in_tree(k,2));
end
addRelativePose(sparse_pg, [1, 0, 0], measInf, 3, 4);

figure
show(sparse_pg,'IDs','off');
title('sparsified pose graph');