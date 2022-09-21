classdef pose2DNode
    %POSE2DNODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pose;           % pose 2D of the node [x, y, theta]
        id;             % ID of the node
        id_graph;       % ID in the pose graph
        lmk_id_list;    % list of the landmarks  
        factors_SE2;    % relative pose factors
        factors_lmk;    % landmark observation factors
    end
    
    methods
        function obj = pose2DNode(pose, id, id_graph)
            %POSE2DNODE Construct an instance of this class
            obj.pose = pose;
            obj.id = id;
            obj.id_graph = id_graph;
            obj.factors_SE2 = {};
            obj.factors_lmk = {};
        end

        function obj = matchLandmark(obj, lmk_id, lmk_factor)
            obj.lmk_id_list(end+1) = lmk_id;
            obj.factors_lmk{end+1} = lmk_factor;
        end

        function obj = addSE2Factor(obj, SE2_factor)
            obj.factors_SE2{end+1} = SE2_factor;
        end

    end
end

