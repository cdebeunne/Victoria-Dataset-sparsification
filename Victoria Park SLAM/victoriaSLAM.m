%%  Load Victoria Dataset data

load("victoriaParkDataset.mat","controllerInput", ...
     "measurements","gpsLatLong","deadReckoning");

poses2D = {};
lmks = {};

%%  Set up parameters
initial_state = [gpsLatLong(1,2) gpsLatLong(1,1) deg2rad(37)]';
initial_covar = eps*eye(3);

% proprioceptive parameters
sigma_velocity = 2;      % [m/s]
sigma_steer= deg2rad(15); % [rad]
process_noise = [sigma_velocity^2 0; 0 sigma_steer^2];
L = 2.83;                % [m]

% exteroceptive parameters
sigma_range = 1;            % [m]
sigma_bearing = deg2rad(3); % [rad]
meas_covar = [sigma_range^2 0; 0 sigma_bearing^2];

max_sensor_range = 30; % [m]

timestep = 0.025; % [sec]

landmark_rejection_thres = 3;  % maximum distance for association
landmark_augmentation_thres = 12; % minimum distance for creation of new landmark
validation_gate = [landmark_rejection_thres landmark_augmentation_thres];

% optimizer parameters 
optimize_step = 2; % Number of update step between each optimization

%%  Loop

victoria_pg = poseGraph;
current_pose = [0 0 0];
last_pose = current_pose;
M = eye(3);

% Update poses cell
count_pose = 1;
count_lmk = 0;
poses2D{count_pose} = pose2DNode(current_pose, count_pose, 1);


for count = 1 : 2500
    % Handle controller input
    controller_input = controllerInput(count, :);
    
    % turn into relative pose
    [relative_pose, info_mat] = carInputFactor(controller_input, L, timestep, process_noise);

    % update odometry integration
    M = M * poseToSE2(relative_pose);
    new_pose = SE2ToPose(poseToSE2(last_pose) * M);
    current_pose = new_pose;
    
    % Handle landmarks
    observed_landmarks = measurements{count};
    if ~isempty(observed_landmarks)
        
        % get current pose
        current_pose = SE2ToPose( poseToSE2(last_pose) * M )';

        disp('----')
        disp(current_pose)

        % Update graph & pose cell
        addRelativePose(victoria_pg, SE2ToPose(M), mat_to_vec(info_mat));
        count_pose = count_pose + 1;
        current_node = victoria_pg.NumNodes;

       
        % Add SE2 factor to previous pose
        poses2D{count_pose} = pose2DNode(current_pose, count_pose, current_node);
        SE2_factor = FactorSE2(SE2ToPose(M), info_mat,...
            count_pose - 1 , count_pose);
        poses2D{count_pose} = addSE2Factor(poses2D{count_pose}, SE2_factor);
        
        
        % Reset odometry integration
        theta = current_pose(3);
        R = [cos(theta), -sin(theta); sin(theta) cos(theta)];
        x = [current_pose(1); current_pose(2)];
        M = eye(3);
        
        % Data association
        for i = 1 : size(observed_landmarks, 1) 
            [l_obs, I_obs] = rangeAndBearing2DFactor(observed_landmarks(i, :), meas_covar);
            I_obs_vec = mat_to_vec(I_obs);
            min_id = 0;
            min_dist = 10000;
            
            % Data association with ML approach
            for j = 1 : count_lmk

                l_known = lmks{j}.coordinates;

                z = R' * (l_known - x);

                if (norm(z) > max_sensor_range)
                    continue;
                end

                mahal_dist = sqrt((l_obs - z)' * I_obs * (l_obs - z));
                
                if (mahal_dist < min_dist)
                    min_dist = mahal_dist;
                    min_id = j;
                end
            end

            % Is it a new landmark ?
            if (min_dist > landmark_augmentation_thres)
                addPointLandmark(victoria_pg, l_obs, I_obs_vec, current_node);
                count_lmk = count_lmk + 1;
                l_world = poseToSE2(current_pose) * [l_obs; 1];
                lmks{count_lmk} = landmarkNode(l_world(1:2), count_lmk, victoria_pg.LandmarkNodeIDs(end));
                
                % add a match
                lmk_factor = FactorLmk(l_obs, I_obs, count_pose, count_lmk);       
                poses2D{count_pose} = matchLandmark(poses2D{count_pose}, count_lmk,lmk_factor);

            % Is it a match ? 
            elseif (min_dist < landmark_rejection_thres)
                addPointLandmark(victoria_pg, l_obs, I_obs_vec, current_node, lmks{min_id}.id_graph);

                % add match
                lmk_factor = FactorLmk(l_obs, I_obs, count_pose, count_lmk);
                poses2D{count_pose} = matchLandmark(poses2D{count_pose}, min_id, lmk_factor);

            end

        end
        if mod(count, optimize_step) == 0
            disp('optim')

            victoria_pg = optimizePoseGraph(victoria_pg, 'builtin-trust-region',...
        'VerboseOutput', 'off',...
        'InitialTrustRegionRadius',10);
            poses2D = updatePose(poses2D, victoria_pg);
            lmks = updateLandmark(lmks, victoria_pg);
        end
        last_pose = poses2D{count_pose}.pose;
        current_pose = last_pose;
    end
end
figure
show(victoria_pg,'IDs','off');
title('victoria graph');
