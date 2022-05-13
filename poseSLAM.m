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

node_indices = (1:updatedPG.NumNodes);
nodes_to_remove = [];
nodes_to_keep = [];

% node pairs that close loopsss
lc_node_pairs = updatedPG.edgeNodePairs(updatedPG.LoopClosureEdgeIDs);

p = 0.5; % probability of removing a node
for i=1:updatedPG.NumNodes

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

%% Full Marginalization

I = computeMatInfJac(updatedPG, node_indices);

% node indice -> place in matrix
map_node = containers.Map(node_indices, node_indices);

% place in matrix -> node indice
map_index = containers.Map(node_indices, node_indices);

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
    map_index(k) = i;
    map_index(i) = k;
    map_node(i) = k;
    map_node(k) = i;
end

% Compute the marginalized Information matrix with Schur Complement
I_aa = I(1:3*n_keep, 1:3*n_keep);
I_bb = I(3*n_keep+1:3*(n_keep+n_remove), 3*n_keep+1:3*(n_keep+n_remove));
I_ab = I(1:3*n_keep,3*n_keep+1:3*(n_keep+n_remove));

I_marg = I_aa - I_ab * inv(I_bb) * I_ab';

% plot it
image = I_marg > 0;
imagesc(image);

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
image = I_marg > 0;
imagesc(image);


