classdef FactorSE2
    %FACTORSE2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        measurement;
        info;
        id1;
        id2;

    end
    
    methods
        function obj = FactorSE2(measurement, info, id1, id2)
            %FACTORSE2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.measurement = measurement;
            obj.info = info;
            obj.id1 = id1;
            obj.id2 = id2;
        end
    end
end

