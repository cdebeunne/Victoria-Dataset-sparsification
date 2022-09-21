function poses2D = updatePose(poses2D, pg)
%UPDATEPOSE Summary of this function goes here
%   Detailed explanation goes here
lmk_graph_ids = pg.LandmarkNodeIDs;

for count_node = 1:pg.NumNodes

    % continue if it is a lmk
    if ismember(count_node, lmk_graph_ids)
        continue
    end

    for count_pose = 1:length(poses2D)

        if (count_node == poses2D{count_pose}.id_graph)
            poses2D{count_pose}.pose = nodeEstimates(pg, count_node);
        end
    end
end

end

