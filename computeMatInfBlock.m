function I_aa = computeMatInfBlock(pg, node_pairs, map)
%
%   Function that extracts a block information matrix from a pose graph
% given the corresponding node pair and its dimensions

I_aa = zeros(3*length(map), 3*length(map));
for k=1:length(node_pairs)
    edge_ij = findEdgeID(pg, node_pairs(k,:));
    [~, Iij_vec] = edgeConstraints(pg, edge_ij);
    Iij = vec_to_mat(Iij_vec);

    i = node_pairs(k,1);
    j = node_pairs(k,2);
    I_aa((map(i)-1)*3+1:map(i)*3,...
        (map(j)-1)*3+1:map(j)*3) = Iij;
    I_aa((map(j)-1)*3+1:map(j)*3,...
        (map(i)-1)*3+1:map(i)*3) = Iij';
end
end

