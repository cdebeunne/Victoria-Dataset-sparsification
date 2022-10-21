classdef VictoriaSLAMInput
    % A class to save an input of the SLAM and to write it in a txt file
    
    properties
        idx;
        type;
        id;
        measurement;
        info_vec;
    end
    
    methods
        function obj = VictoriaSLAMInput(timestamp, type, id, measurement, info_mat)
            % Constructor
            obj.idx = timestamp;
            obj.type = type;
            obj.id = id;
            obj.measurement = measurement;
            obj.info_vec = mat_to_vec(info_mat);
        end

        function str = writeInput(obj)
            % The line to write in a text file
            if (obj.type == "landmark")
                str = (obj.idx + " , " + obj.type + " , " +  obj.id + " , " + ...
                    obj.measurement(1) + " , " + obj.measurement(2) + " , "  +...
                    obj.info_vec(1) + " , " + obj.info_vec(2) + " , " + obj.info_vec(3));
            elseif (obj.type == "odometry")
                str = (obj.idx + " , " + obj.type + " , " +  ...
                    obj.measurement(1) + " , " + obj.measurement(2) + " , "  + obj.measurement(3) + " , " + ...
                    obj.info_vec(1) + " , " + obj.info_vec(4) + " , " + obj.info_vec(6));
            end
        end
  
    end
end

