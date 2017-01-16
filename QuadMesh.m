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
                dof=reshape([obj.NN(obj.NEL(i,:),4),obj.NN(obj.NEL(i,:),5)]',[8,1]);
                obj.Elements(i)=Element(x,y,gNodes,dof,obj.orderInt);
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
        
        function actualizedPoints=createPoints(obj,numPoints,insertedPoints)
            % Creates a Random grid of points starting at the the origin
            % and right and up from there.
            % insertedPoints is a list of Points that are inserted
            meshLimits=obj.meshLimits;
            actualizedPoints=insertedPoints;
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
        
        function WeightatQuadrature(obj,Cloud)
           % Function takes in cloud handle to calculate weight functions 
           size1=Cloud.numberOfNodes;
           size2=obj.orderInt^2;
           Elements=obj.Elements; % Create Module Copy of Element Data Structure
           parfor i=1:obj.noElements
               wtemp=zeros(3,size1,size2);
               myElement=Elements(i); % Create Loop Copy of Array
               for j=1:size2
                   Cords=myElement.getIntCord(j); x=Cords(1); y=Cords(2);
                    for a=1:size1
                       CordsA=Cloud.Nodes(a).cordinates; 
                       if (CordsA(1)-x)<=Cloud.Nodes(a).a && (CordsA(2)-y)<=Cloud.Nodes(a).a
                           weightstemp=Cloud.Nodes(a).weight.weights([x;y]);
                           wtemp(1,a,j)=Cloud.Nodes(a).weight.w([x;y],weightstemp);
                           wxtemp=Cloud.Nodes(a).weight.wx([x;y],weightstemp);
                           wtemp(2,a,j)=wxtemp(1);
                           wtemp(3,a,j)=wxtemp(2);
                       end
                    end
               end
               myElement.setStoredWeight(wtemp);
               Elements(i)=myElement; % Assign Local Handle to Module
           end
           obj.Elements=Elements; % Assign Module to Global
        end
        
        function [RKv,RKdx,RKdy]=RKPMatQuadrature(obj,Cloud)
           % Function takes in Cloud handle to calculate RKPM shape
           % function values at each Quadrature Point for later integration
           size1=Cloud.numberOfNodes;
           size2=obj.orderInt^2;
           RKv=zeros(size1,obj.noElements,size2,1);
           RKdx=RKv;
           RKdy=RKv;
           Elements=obj.Elements; % Create Module Copy of Element Data Structure
           for i=1:obj.noElements
               myElement=Elements(i); % Create Loop Copy of Array
               StoredRKPM=zeros(3,size1,size2);
               weightsAtQuadrature=myElement.StoredWeight;
               for j=1:size2
                   Cords=myElement.getIntCord(j); x=Cords(1); y=Cords(2);
                   wAQ=weightsAtQuadrature(:,:,j);
                   for a=1:size1
                       %CordsA=Cloud.Nodes(a).cordinates; 
                       if wAQ(1,a)~=0
                           StoredRKPM(1,a,j)=Cloud.Nodes(a).sF.getValue([x;y],wAQ,a);
                           DX=Cloud.Nodes(a).sF.getValueDx([x;y],wAQ,a);
                           StoredRKPM(2,a,j)=DX(1);
                           StoredRKPM(3,a,j)=DX(2);
                           RKv(a,i,j,1)=StoredRKPM(1,a,j);
                           RKdx(a,i,j,1)=DX(1);
                           RKdy(a,i,j,1)=DX(2);
                       else
                           %Pass
                       end
                   end
                   myElement.setStoredRKPM(StoredRKPM);           
                   Elements(i)=myElement; % Assign Local Handle to Module
               end
           end
           obj.Elements=Elements; % Assign Module to Global
        end
        
    end
end

