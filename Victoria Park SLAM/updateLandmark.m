function lmks = updateLandmark(lmks, pg)
%UDPATELANDMARK Summary of this function goes here
%   Detailed explanation goes here
lmk_graph_ids = pg.LandmarkNodeIDs;

for count_node = 1:length(lmk_graph_ids)
    id = lmk_graph_ids(count_node);

    for count_lmk = 1:length(lmks)

        if (id == lmks{count_lmk}.id_graph)
            coord = nodeEstimates(pg, id);
            lmks{count_lmk}.coordinates = coord(1:2)';
        end
    end
end
end

