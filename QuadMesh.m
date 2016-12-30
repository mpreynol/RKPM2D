classdef QuadMesh < handle
    %QUADMESH is a quadrature mesh for RKPM in 2D
    
    properties 
        orderInt; % level of quadrature for mesh
        NN; % Array of Grid Points
        NEL; % Element Connectivity
        noElements; % No. for elements
        Elements; % Collection of Element Objects
    end
    
    methods
        function obj=QuadMesh(NN,NEL,orderInt)
            obj.orderInt=orderInt;
            obj.NN=NN;
            obj.NEL=NEL;
            obj.noElements=size(NEL,1);
            obj.Elements=Element.empty(obj.noElements,0); % Declare Empty Collection
            obj.setElements()
            
        end
        function setElements(obj)
            for i=1:obj.noElements
                gNodes=obj.NEL(i,:); x=obj.NN(gNodes,2); y=obj.NN(gNodes,3);
                obj.Elements(i)=Element(x,y,gNodes,obj.orderInt);
            end
        end
    end
    
end

