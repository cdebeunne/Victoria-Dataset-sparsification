function pg = load_g2o(g2o_data_file)
%
% function measurements = load_g2o(g2o_data_file)
% This function accepts as input a .g2o file containing the description of
% a 2D or 3D pose graph SLAM problem, and returns a MATLAB poseGraph

fid = fopen(g2o_data_file, 'r');
edge_id = 0;
read_line = fgets(fid);  % Read the next line from the file
pg = poseGraph;


while ischar(read_line)  % As long as this line is a valid character string
    
    token = strtok(read_line);

    if(strcmp(token, 'EDGE_SE2'))
        % 2D OBSERVATION
        
        edge_id = edge_id + 1;
        
        % The g2o format specifies a 3D relative pose measurement in the
        % following form:
        
        % EDGE_SE2 id1 id2 dx dy dtheta, I11, I12, I13, I22, I23, I33
        
        [~, id1, id2, dx, dy, dth, I11, I12, I13, I22, I23, I33] = ...
            strread(read_line, '%s %d %d %f %f %f %f %f %f %f %f %f');
        
         % Store the connectivity of this edge
        connectivity = [id1 + 1, id2 + 1];  % NB: .g2o uses 0-based indexing, whereas MATLAB uses 1-based indexing
        
        % Store the translational measurement
        measurement = [dx, dy, dth]';
        
        % Reconstruct the information matrix
        measurement_info = [I11, I12, I13, I22, I23, I33];

        addRelativePose(pg, measurement, measurement_info, connectivity(1), connectivity(2));
        
       
    end
    
    read_line = fgets(fid);
end


end

