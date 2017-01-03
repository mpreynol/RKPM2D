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
                obj.Nodes(i)=RKNode(obj.nodalData(i,1),[obj.nodalData(i,2);obj.nodalData(i,3)],obj.nodalData(i,4),obj.order,obj,obj.nodalData(i,5));
            end
        end
        
        function A=distToAll(obj)
           % Matrix with each element a data structure with distance vector between points i and j (sym).
           A=zeros(obj.numberOfNodes,obj.numberOfNodes,3);
           for i=1:obj.numberOfNodes
              for j=1:obj.numberOfNodes
                  if i~=j
                        A(i,j,1)=obj.nodalData(i,2)-obj.nodalData(j,2);
                        A(i,j,2)=obj.nodalData(i,3)-obj.nodalData(j,3);
                        A(i,j,3)=sqrt(A(i,j,1)^2+A(i,j,2)^2);
                        A(i,j,1)=A(i,j,1)./ A(i,j,3);
                        A(i,j,2)=A(i,j,2)./ A(i,j,3);
                   end
              end
           end
        end
    end
    
end

