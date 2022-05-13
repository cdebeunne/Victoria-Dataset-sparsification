%% Create a toy example

% We create a toy pose graph to understand basic concept
%
%                      > x4 ---> x5
%                    /    ^
%                   /     |  
%       x1 ---> x2 --->  x3

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
for k = 1:n_remove
    markov_blanket_remove = [markov_blanket_remove,...
                           getMarkovBlanket(updated_toy_pg, nodes_to_remove(k))];
end
markov_blanket_remove = unique(markov_blanket_remove);
n_mb = length(markov_blanket_remove);

% here the ordering is the same as markov_blanket_keep
I_t = computeMatInfJac(updated_toy_pg, markov_blanket_remove);

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
image = I_marg > 0;
imagesc(image);

