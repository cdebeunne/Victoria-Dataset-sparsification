classdef FactorLmk
    %FACTORLMK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        measurement;
        info;
        id_pose;
        id_lmk;
    end
    
    methods
        function obj = FactorLmk(measurement, info, id_pose, id_lmk)
            %FACTORLMK Construct an instance of this class
            %   Detailed explanation goes here
            obj.measurement = measurement;
            obj.info = info;
            obj.id_pose = id_pose;
            obj.id_lmk = id_lmk;
        end
    end
end

