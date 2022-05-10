function markov_blanket = getMarkovBlanket(pg,node)
% Returns the nodes in the markov blanket of a node in a given pose graph

node_pairs = pg.edgeNodePairs;
markov_blanket = [];

for k = 1:length(node_pairs)
    if (node_pairs(k,1) == node)
        markov_blanket = [markov_blanket, node_pairs(k,2)];
    elseif (node_pairs(k,2) == node)
        markov_blanket = [markov_blanket, node_pairs(k,1)];
    end
end


end

