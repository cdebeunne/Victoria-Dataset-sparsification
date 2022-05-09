function Lambda = computeMatInf(pg)
% This function returns the information matrix of a SE2 posegraph

Lambda = zeros(pg.NumNodes*3, pg.NumNodes*3);
nodePairs = edgeNodePairs(pg);
for i= 1:length(nodePairs)
    edgeij = findEdgeID(pg,nodePairs(i,:));
    [~, Lambdaij_vec] = edgeConstraints(pg, edgeij);
    Lambdaij = [Lambdaij_vec(1,1), Lambdaij_vec(1,2), Lambdaij_vec(1,3);
        Lambdaij_vec(1,2), Lambdaij_vec(1,4), Lambdaij_vec(1,5);
        Lambdaij_vec(1,3), Lambdaij_vec(1,5), Lambdaij_vec(1,6)];
    Lambda((nodePairs(i,1)-1)*3+1:nodePairs(i,1)*3, (nodePairs(i,2)-1)*3+1:nodePairs(i,2)*3) = Lambdaij;
    Lambda((nodePairs(i,2)-1)*3+1:nodePairs(i,2)*3, (nodePairs(i,1)-1)*3+1:nodePairs(i,1)*3) = Lambdaij';
end
end

