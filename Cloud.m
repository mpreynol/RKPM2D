classdef Cloud < handle
    %CLOUD is a list of node objects for an RKPM solution space
    % Currently using 2D
    
    properties
        nodalData=[]; % Array of nodal Cordinates
        Nodes=[]; % List of Object Nodes
        numberOfNodes; % Number of nodes in the cloud
        order; % Order of basis Functions
    end
    
    methods
        function obj=Cloud(nodalData,order)
            %Constructor
            obj.nodalData=nodalData;
            obj.numberOfNodes=size(nodalData,1);
            obj.order=order;
            obj.setNodeObjects();
            
        end
        
        function setNodeObjects(obj)
            obj.Nodes=RKNode.empty(obj.numberOfNodes,0);
            for i=1:obj.numberOfNodes
                obj.Nodes(i)=RKNode(obj.nodalData(i,1),[obj.nodalData(i,2),obj.nodalData(i,3)],obj.nodalData(i,4),obj.order,obj,obj.nodalData(i,5));
            end
        end
    end
    
end

