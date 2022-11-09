function Lambda = computeMatInfJac(pg, nodes)
%   Function that extracts a block information matrix from a pose graph
% given the set of nodes of interest 

node_pairs = edgeNodePairs(pg);
n_lmk = sum(ismember(pg.LandmarkNodeIDs, nodes));
n_pose = length(nodes) - n_lmk;
dim_lambda = n_pose*3 + n_lmk*2;

Lambda = zeros(dim_lambda, dim_lambda);

for k = 1:length(node_pairs)

    if ~(ismember(node_pairs(k,1), nodes) && ...
         ismember(node_pairs(k,2), nodes))
       continue 
    end
    
    % Retrieve the factor id and information matrix
    edge_ij = findEdgeID(pg, node_pairs(k,:));
    [~, Iij_vec] = edgeConstraints(pg, edge_ij);
    Iij = vec_to_mat(Iij_vec);
    i = node_pairs(k,1);
    j = node_pairs(k,2);
    
    % Case of relative pose factor
    if (size(Iij, 1) == 3)

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
        H = zeros(3, dim_lambda);
        H(:, (p-1)*3+1:p*3) = -[Rj'*Ri Rperp*Rj'*(ti-tj);
                                        0 0 1];
        H(:, (q-1)*3+1:q*3) = eye(3,3);
    
        Lambda = Lambda + H'*Iij*H;
    end 

    % Case of a landmark factor 
    if (size(Iij, 1) == 2)

        % retrieve the robot pose
        meas_i = nodeEstimates(pg,i);
        s = sin(meas_i(3));
        c = cos(meas_i(3));
        x = meas_i(1:2);
        R = [c -s;
            s c];

        % retrieve the landmark pose
        meas_j = nodeEstimates(pg,j);
        l = meas_j(1:2);


        % Retrieve indices in the information matrix
        p = find(nodes==i);
        q = find(nodes==j);

        % s_q is the index in the information matrix of the landmark
        s_q = n_pose * 3 + (q-n_pose - 1) * 2 + 1;

        % Compute 2D jacobians
        % Formula from debeunne et al that is probably wrong
        H = zeros(2, dim_lambda);
        dz_dtheta = [s * x(1) - c * x(2); c * x(1) + s * x(2)] +...
            [-s * l(1) + c * l(2); -c * l(1) - s * l(2)];
        dz_dx = -R';
        H(:, (p-1)*3+1:p*3) = [dz_dx dz_dtheta];
        H(:, s_q:s_q+1) = R';
        Lambda = Lambda + H'*Iij*H;
    end
end
end

