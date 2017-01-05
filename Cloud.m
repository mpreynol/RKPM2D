classdef Cloud < handle
    %CLOUD is a list of node objects for an RKPM solution space
    % Currently using 2D
    
    properties
        nodalData=[]; % Array of nodal Cordinates
        Nodes=[]; % List of Object Nodes
        numberOfNodes; % Number of nodes in the cloud
        order; % Order of basis Functions
        u; % Solutions *** Lack of Kronecker Means this isn't the PDE Solution Yet
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
                obj.Nodes(i)=RKNode(obj.nodalData(i,1),[obj.nodalData(i,2);obj.nodalData(i,3)],obj.nodalData(i,4),obj.order,[obj.nodalData(i,5);obj.nodalData(i,6)],obj);
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
        
        function plotCloud(obj)
            for i=1:obj.numberOfNodes
                refresh
                scatter(obj.nodalData(i,2),obj.nodalData(i,3),'ko')
                hold on
            end
        end
        
        function [ufull] = reAssembleUnknowns(obj,ureduced,BE)
            % Method reassembles a full 'u' vector for a reduced oned
            L=BE==-inf;
            ufull=zeros(length(BE),1);
            counter=1;
            for i=1:length(BE)
               if L(i)==1
                   ufull(i)=ureduced(counter);
                   counter=counter+1;
               else
                   ufull(i)=BE(i);
               end
            end
        end
        
        function parseSolution(obj,u)
            % Method insets solved for coeff's 'u' into the Cloud and Nodes
            obj.u=u;
            for i=1:obj.numberOfNodes
                dof=[obj.Nodes(i).nodeNumber*2-1;obj.Nodes(i).nodeNumber*2];
                obj.Nodes(i).u=u(dof);
            end
        end
        
        function plotCloudU(obj,scale)
                for j=1:obj.numberOfNodes
                    Cords=obj.Nodes(j).cordinates;
                    deformation=0;
                    for i=1:obj.numberOfNodes
                        deformation=deformation+obj.Nodes(i).sF.getValue([Cords(1);Cords(2)])*obj.Nodes(i).u;
                    end
                    scatter(Cords(1)+scale*deformation(1),Cords(2)+scale*deformation(2),'bo')
                    hold on
                end
        end
        
        function [K,F]=integrateDomain(obj,C,Q,Mesh)
            n=obj.numberOfNodes;
            K=zeros(n*2);
            F=zeros(n*2,1);
            totalPoints=Mesh.noElements*size(Mesh.Elements(1).G2,1);
            quadCounter=1;
            for k=1:Mesh.noElements % Loop through Elements
                for l=1:size(Mesh.Elements(k).G2,1) % Loop through Quadrature Points
                    display([num2str(quadCounter * 100 / totalPoints) , ' Percent Complete'])
                    quadCounter=quadCounter+1;
                    Cords=Mesh.Elements(k).getIntCord(l); x=Cords(1); y=Cords(2); Jw=Cords(3); % Jw = Jacobian * Quadrature Weight
                    for a=1:n % Loop over Shape Functions
                        for b=1:n % Loop over Shape Functions
                            if b>=a
                                CordsA=obj.Nodes(a).cordinates; CordsB=obj.Nodes(b).cordinates;
                                if ((CordsA(1)-CordsB(1))^2+(CordsA(2)-CordsB(2))^2)<obj.Nodes(a).a^2
                                    NDa=obj.Nodes(a).sF.getValueDx([x;y]);
                                    NDb=obj.Nodes(b).sF.getValueDx([x;y]);
                                    Ba=[NDa(1),0;0,NDa(2);NDa(2),NDa(1)];
                                    Bb=[NDb(1),0;0,NDb(2);NDb(2),NDb(1)];
                                    Kab=Ba'*C*Bb*Jw;
                                    K(2*a-1,2*b-1)=K(2*a-1,2*b-1)+Kab(1,1);
                                    K(2*a-1,2*b)=K(2*a-1,2*b)+Kab(1,2);
                                    K(2*a,2*b-1)=K(2*a,2*b-1)+Kab(2,1);
                                    K(2*a,2*b)=K(2*a,2*b)+Kab(2,2);
                                else
                                    %Pass
                                end
                            end
                        end
                        F(2*a-1)=F(2*a-1)+obj.Nodes(a).sF.getValue([x;y])*Jw*Q(1);
                        F(2*a)=F(2*a)+obj.Nodes(a).sF.getValue([x;y])*Jw*Q(2);
                    end
                end
            end
            K=triu(K)+tril(K',-1);
        end
        
        function U=returnInterpolatedU(obj,x)
            U=zeros(2,1);
           for i=1:obj.numberOfNodes
               U=U+obj.Nodes(i).sF.getValue(x)*obj.Nodes(i).u;
           end
        end
    end
    
end

