% Load posegraph
pg = load_g2o('input_MITb_g2o.g2o');
disp(pg);
show(pg,'IDs','off');
title('Original Pose Graph');

% Solve posegraph
updatedPG = optimizePoseGraph(pg);
figure
show(updatedPG,'IDs','off');
title('Updated Pose Graph');

%% Select nodes to keep and to remove 

node_indices = (1:updatedPG.NumNodes);
nodes_to_remove = [];
nodes_to_keep = [];

% node pairs that close loopsss
lc_node_pairs = updatedPG.edgeNodePairs(updatedPG.LoopClosureEdgeIDs);

p = 0.3; % probability of removing a node
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

I = computeMatInfJac(updatedPG);

% node indice -> place in matrix
map_node = containers.Map(node_indices, node_indices);

% place in matrix -> node indice
map_index = containers.Map(node_indices, node_indices);

for k = 1:length(nodes_to_keep)
    % I copy ?
    i = map_node(nodes_to_keep(k));
    I((k-1)*3+1:k*3, :) = I((i-1)*3+1:i*3, :);
    I(:,(k-1)*3+1:k*3) = I(:, (i-1)*3+1:i*3);
    map_node(k) = i;
    map_node(i) = k;
end

% Compute the marginalized Information matrix with Schur Complement
I_aa = I(1:3*n_keep, 1:3*n_keep);
I_bb = I(3*n_keep+1:3*(n_keep+n_remove), 3*n_keep+1:3*(n_keep+n_remove));
I_ab = I(1:3*n_keep,3*n_keep+1:3*(n_keep+n_remove));

I_marg = I_aa - I_ab * inv(I_bb) * I_ab';

% plot it
image = I_marg > 0;
imagesc(image);

%% Compute KLD divergence between pose graphs

nodePairs = edgeNodePairs(pg);
KLD = 0;
for i=1:length(nodePairs)
    edgeij = findEdgeID(pg,nodePairs(i,:));
    [mu, Lambda] = edgeConstraints(pg, edgeij);
    Lambda = vecToMat(Lambda);
    [mu_tilde, Lambda_tilde] = edgeConstraints(updatedPG, edgeij);
    Lambda_tilde = vecToMat(Lambda_tilde);
    KLD = KLD + 0.5*(trace(Lambda_tilde*inv(Lambda))...
                    - log(norm(Lambda_tilde*inv(Lambda)))...
                    + (mu - mu_tilde)*inv(Lambda_tilde)*(mu - mu_tilde)' ...
                    - length(mu));
end
