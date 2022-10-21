load('VictoriaDatasetxylmk.mat')

% The pose graph object from matlab
victoria_pg = poseGraph;

% A container to link the landmarks ids to their ids in the graph
map_lmkid_graphid = containers.Map();

% Loop through each measurements
for i = 1:length(SLAMInputs)
    measurement = SLAMInputs{i};

    if (measurement.type == "odometry")
        addRelativePose(victoria_pg, measurement.measurement', measurement.info_vec);
        current_node = victoria_pg.NumNodes;
    end

    if (measurement.type == "landmark")
        if (isKey(map_lmkid_graphid, num2str(measurement.id)))
            id_lmk = map_lmkid_graphid(num2str(measurement.id));
            addPointLandmark(victoria_pg, measurement.measurement', measurement.info_vec, current_node, id_lmk);
        else
            addPointLandmark(victoria_pg, measurement.measurement', measurement.info_vec);
            map_lmkid_graphid(num2str(measurement.id)) = victoria_pg.LandmarkNodeIDs(end);
        end
    end
end

% Graph optimization
victoria_pg = optimizePoseGraph(victoria_pg, 'builtin-trust-region',...
        'VerboseOutput', 'off',...
        'InitialTrustRegionRadius',10);

figure
show(victoria_pg,'IDs','off');
title('victoria graph');