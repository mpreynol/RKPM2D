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
            % Constructor
            obj.orderInt=orderInt;
            obj.NN=NN;
            obj.NEL=NEL;
            obj.noElements=size(NEL,1);
            obj.Elements=Element.empty(obj.noElements,0); % Declare Empty Collection
            obj.setElements()
            
        end
        function setElements(obj)
            % Loop through and define elements
            for i=1:obj.noElements
                gNodes=obj.NEL(i,:); x=obj.NN(gNodes,2); y=obj.NN(gNodes,3);
                obj.Elements(i)=Element(x,y,gNodes,obj.orderInt);
            end
        end
        
        function ret=insideMesh(obj,x)
            % Function returns if point 'x' is inside Element mesh or not.
            ret=false;
            for i=1:obj.noElements       
                if (sum(x(1)<=obj.Elements(i).x)>0 && sum(x(1)>=obj.Elements(i).x)>0 && ...
                   sum(x(2)<=obj.Elements(i).y)>0 && sum(x(2)>=obj.Elements(i).y)>0) 
                    ret=true;
                end
            end
        end
        
        function area=meshArea(obj)
           % Returns the area of the mesh 
           area=0;
           for i=1:obj.noElements
               for j=1:size(obj.Elements(i).G2)
                   intPoint=[obj.Elements(i).G2(j,1);obj.Elements(i).G2(j,1)];
                   obj.Elements(i).Shape.setAll(intPoint(1),intPoint(2));
                   area=area+obj.Elements(i).Shape.j*obj.Elements(i).G2(j,3);
               end
           end
        end
        
        function meshLimits=meshLimits(obj)
            % Returns the overall mesh limits
            xmin=obj.Elements(1).x(1);
            xmax=obj.Elements(1).x(2);
            ymin=obj.Elements(1).y(1);
            ymax=obj.Elements(1).y(3);
            for i=1:obj.noElements
                if min(obj.Elements(i).x)<xmin
                    xmin=min(obj.Elements(i).x);
                end
                if max(obj.Elements(i).x)>xmax
                    xmax=max(obj.Elements(i).x);
                end
                if min(obj.Elements(i).y)<ymin
                    ymin=min(obj.Elements(i).y);
                end
                if max(obj.Elements(i).y)>ymax
                    ymax=max(obj.Elements(i).y);
                end
            end
            meshLimits=[xmin,xmax,ymin,ymax];
        end
        
        function actualizedPoints=createPoints(obj,numPoints)
            % Creates a Random grid of points starting at the the origin
            % and right and up from there.
            meshLimits=obj.meshLimits;
            actualizedPoints=[];
            Area=obj.meshArea();
            IdealPointDensity=sqrt(Area/(numPoints));
            while size(actualizedPoints,1)<numPoints
                % Add Point
                tempPoint=[rand()*meshLimits(2),rand()*meshLimits(4)];
                % Check Point Location
                minDist=inf;
                for i=1:size(actualizedPoints,1)
                    tempDist=norm(tempPoint-actualizedPoints(i,:));
                    if tempDist<minDist
                        minDist=tempDist;
                    end
                end
                
                BorderTolerance=IdealPointDensity*0.5;
                if abs(tempPoint(1)-meshLimits(1))<BorderTolerance
                    tempPoint(1)=meshLimits(1);
                end
                if abs(tempPoint(1)-meshLimits(2))<BorderTolerance
                    tempPoint(1)=meshLimits(2);
                end
                if abs(tempPoint(2)-meshLimits(3))<BorderTolerance
                    tempPoint(2)=meshLimits(3);
                end
                if abs(tempPoint(2)-meshLimits(4))<BorderTolerance
                    tempPoint(2)=meshLimits(4);
                end
                if minDist>=IdealPointDensity*0.8
                    actualizedPoints=[actualizedPoints;tempPoint];
                end
            end
        end
        
        function RKPMatQuadrature(obj,Cloud)
           % Function takes in Cloud handle to calculate RKPM shape
           % function values at each Quadrature Point
           totalPoints=obj.noElements*obj.orderInt^2;
           quadCounter=1;
           for i=1:obj.noElements
               StoredRKPM=zeros(3,Cloud.numberOfNodes,obj.orderInt^2);
               for j=1:obj.orderInt^2
                   Cords=obj.Elements(i).getIntCord(j); x=Cords(1); y=Cords(2);
                   display([num2str(quadCounter * 100 / totalPoints) , ' Percent Complete'])
                   quadCounter=quadCounter+1;
                   for a=1:Cloud.numberOfNodes
                       CordsA=Cloud.Nodes(a).cordinates; 
                       if ((CordsA(1)-x)^2+(CordsA(2)-y)^2)<Cloud.Nodes(a).a^2
                           StoredRKPM(1,a,j)=Cloud.Nodes(a).sF.getValue([x;y]);
                           DX=Cloud.Nodes(a).sF.getValueDx([x;y]);
                           StoredRKPM(2,a,j)=DX(1);
                           StoredRKPM(3,a,j)=DX(2);
                       else
                           %Pass
                       end
                   end
                   obj.Elements(i).setStoredRKPM(StoredRKPM);           
               end
           end
        end
        
    end
end

