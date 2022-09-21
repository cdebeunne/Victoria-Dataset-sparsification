classdef landmarkNode
    %LANDMARKNODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coordinates;        % Cartesian 2D coordinates of the lmk
        id;                 % Id of the lmk
        id_graph;           % Id in the graph
    end
    
    methods
        function obj = landmarkNode(coordinates, id, id_graph)
            %LANDMARKNODE Construct an instance of this class
            %   Detailed explanation goes here
            obj.coordinates = coordinates;
            obj.id = id;
            obj.id_graph = id_graph;
        end
        
    end
end

