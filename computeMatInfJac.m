function Lambda = computeMatInfJac(pg, nodes)
%
%   Function that extracts a block information matrix from a pose graph
% given the set of nodes of interest 

Lambda = zeros(length(nodes)*3, length(nodes)*3);
node_pairs = edgeNodePairs(pg);

for k = 1:length(node_pairs)

    if ~(ismember(node_pairs(k,1), nodes) && ...
         ismember(node_pairs(k,2), nodes))
       continue 
    end

    edge_ij = findEdgeID(pg, node_pairs(k,:));
    [~, Iij_vec] = edgeConstraints(pg, edge_ij);
    Iij = vec_to_mat(Iij_vec);

    i = node_pairs(k,1);
    j = node_pairs(k,2);
    
    % Retrieve node estimates
    meas_i = nodeEstimates(pg,i);
    ti = meas_i(1:2)';
    Ri = [cos(meas_i(3)) -sin(meas_i(3));
        sin(meas_i(3)) cos(meas_i(3))];

    meas_j = nodeEstimates(pg,j);
    tj = meas_j(1:2)';
    Rj = [cos(meas_j(3)) -sin(meas_j(3));
        sin(meas_j(3)) cos(meas_j(3))];
    Rperp = [0 1; -1 0];

    % Retrieve indices in the information matrix
    p = find(nodes==i);
    q = find(nodes==j);

    % Compute 2D jacobians
    % Formula from 2D poseSLAM in GTSAM Dellaert.
    H = zeros(3, 3*length(nodes));
    H(:, (p-1)*3+1:p*3) = [Rj'*Ri Rperp*Rj'*(ti-tj);
                                    0 0 1];
    H(:, (q-1)*3+1:q*3) = eye(3,3);

    Lambda = Lambda + H'*Iij*H;
end
end

