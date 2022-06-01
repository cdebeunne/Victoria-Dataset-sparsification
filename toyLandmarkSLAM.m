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

% TO IMPROVE
% add a map to the type of node to build an information matrix with
% different node sizes
pose_nodes = [1 2 3];
lmk_nodes = [4 5];

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


% select the nodes to keep in the markov blanket
nodes_to_keep_t = nodes_to_keep(ismember(nodes_to_keep,markov_blanket_remove));
n_keep_t = length(nodes_to_keep_t);

% We select the rows and the columns to keep and to remove
rows_to_keep = [];
rows_to_remove = [];
for k = 1:length(nodes_to_keep_t)
    j = nodes_to_keep_t(k);
    rows_to_keep = [rows_to_keep, (j-1)*3+1, (j-1)*3+2,...
        (j-1)*3+3];
end
for k = 1:length(node_to_remove)
    j = node_to_remove(k);
    rows_to_remove = [rows_to_remove, (j-1)*3+1, (j-1)*3+2,...
        (j-1)*3+3];
end

% Compute the marginalized Information matrix with Schur Complement
I_aa = I_t(rows_to_keep, rows_to_keep);
I_bb = I_t(rows_to_remove, rows_to_remove);
I_ab = I_t(rows_to_keep, rows_to_remove);

I_marg = I_aa - I_ab * pinv(I_bb) * I_ab';

% plot it
imagesc(I_marg > 0);
title("Marginalized, Dense Information matrix on elimination clique")

%% Select Manually a topology for the new subgraph

%            > l4 <   
%          /        \  
%         /          \  
%    |--x1 --------> x3
%         \          /
%          \        /
%            > l5 < 
% This is where there is room for research 
% BEWARE the closed form solution cannot be used with loops!

edges_in_tree = [map_node_t(1) map_node_t(4);
                 map_node_t(1) map_node_t(5);
                 map_node_t(3) map_node_t(4);
                 map_node_t(3) map_node_t(5);
                 map_node_t(1) map_node_t(3)];

