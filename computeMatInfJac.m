function Lambda = computeMatInfJac(pg)
%
%   Function that extracts a block information matrix from a pose graph
% given the corresponding node pair and its dimensions

Lambda = zeros(pg.NumNodes*3, pg.NumNodes*3);
node_pairs = edgeNodePairs(pg);
for k = 1:length(node_pairs)
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
    
    % Compute 2D jacobians
    % Formula from 2D poseSLAM in GTSAM Dellaert.
    H = zeros(3, 3*pg.NumNodes);
    H(:, (i-1)*3+1:i*3) = [Rj'*Ri Rperp*Rj'*(ti-tj);
                                    0 0 1];
    H(:, (j-1)*3+1:j*3) = eye(3,3);

    Lambda = Lambda + H'*Iij*H;
end
end

