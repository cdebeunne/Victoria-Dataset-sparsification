%% Create a toy example

% We create a toy pose graph to understand basic concept
%
%            > l4     > l5
%          /    ^   /    ^
%         /     |  /     | 
%    |--x1 ---> x2 ---> x3
%
% The observation model is z = R' * (l - x)

toy_pg = poseGraph;
meas_inf_proprio = [1.5 0 0 4 0 400];
meas_inf_extero = [10 0 10];

% Relative pose
addRelativePose(toy_pg, [1, 0, pi/4], meas_inf_proprio, 1, 2);
addRelativePose(toy_pg, [0.7071, 0, 0], meas_inf_proprio, 2, 3);

% Landmark measurement
addPointLandmark(toy_pg, [0.5 0.5], meas_inf_extero, 1, 4);
addPointLandmark(toy_pg, [-0.5, 0.5], meas_inf_extero, 2, 4);
addPointLandmark(toy_pg, [0.3536, 0.3536], meas_inf_extero, 2, 5);
addPointLandmark(toy_pg, [-0.3536, 0.3536], meas_inf_extero, 3, 5);


figure
show(toy_pg,'IDs','off');
title('toy graph');

%% Choose the node we want to remove and bookeping

node_to_remove = 2;
nodes_to_keep = [1 3 4 5];

%% Marginalization

% get the global markov blanket 
markov_blanket_remove = [];
markov_blanket_remove = [markov_blanket_remove,...
                       getMarkovBlanket(toy_pg, node_to_remove)];
markov_blanket_remove = unique(markov_blanket_remove);
n_mb = length(markov_blanket_remove);

% compute the information matrix on the elimination clique
I_t = computeMatInfJac(toy_pg, markov_blanket_remove);

% we need to add a unary factor so that I_t is full rank
H = zeros(3, 3*n_mb);
H(1:3, 1:3) = eye(3,3);
I = 100 * eye(3,3);
I_t = I_t + H' * I * H;

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
title("Marginalized, Dense Information matrix on elimination clique")