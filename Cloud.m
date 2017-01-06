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
                scatter(obj.nodalData(i,2),obj.nodalData(i,3),'ko','filled')
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
            for k=1:Mesh.noElements % Loop through Elements
                for l=1:size(Mesh.Elements(k).G2,1) % Loop through Quadrature Points
                    Cords=Mesh.Elements(k).getIntCord(l); Jw=Cords(3); % Jw = Jacobian * Quadrature Weight
                    for a=1:n % Loop over Shape Functions
                        for b=1:n % Loop over Shape Functions
                            if b>=a
                                CordsA=obj.Nodes(a).cordinates; CordsB=obj.Nodes(b).cordinates;
                                if ((CordsA(1)-CordsB(1))^2+(CordsA(2)-CordsB(2))^2)<obj.Nodes(a).a^2
                                    NDa=Mesh.Elements(k).StoredRKPM(2:3,a,l);
                                    NDb=Mesh.Elements(k).StoredRKPM(2:3,b,l);
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
                        if sum(Q)>0
                            F(2*a-1)=F(2*a-1)+Mesh.Elements(k).StoredRKPM(1,a,l)*Jw*Q(1);
                            F(2*a)=F(2*a)+Mesh.Elements(k).StoredRKPM(1,a,l)*Jw*Q(2);
                        end
                    end
                end
            end
            K=triu(K)+tril(K',-1);
        end
        
        function F=integrateBoundary(obj,Mesh,BN)
            F=zeros(2*obj.numberOfNodes,1);
            for j=1:Mesh.noElements
                h=BN(Mesh.Elements(j).dof);              
                if sum((sum(h~=0))) % Then we have tractions on te element
                    nInt=Mesh.Elements(j).orderInt;
                    G1=Mesh.Elements(j).G1;
                    if (h(1)~=0 && h(3)~=0) || (h(2)~=0 && h(4)~=0)
                        % Surface 1: eta=-1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(G1(i,1),-1);
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xxi)^2+(Mesh.Elements(j).Shape.Yxi)^2);
                            hInt=[h(1);h(2)]*(1-(1+G1(i,1))/2)+(1+G1(i,1))/2*[h(3);h(4)];
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*hInt*js*G1(i,2);
                                end
                            end
                        end
                    end
                    if (h(3)~=0 && h(5)~=0) || (h(4)~=0 && h(6)~=0)
                        % Surface 2: xi=1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(1,G1(i,1));
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xeta)^2+(Mesh.Elements(j).Shape.Yeta)^2);
                            hInt=[h(3);h(4)]*(1-(1+G1(i,1))/2)+(1+G1(i,1))/2*[h(5);h(6)];
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*hInt*js*G1(i,2);
                                end
                            end
                        end
                    end
                    if (h(5)~=0 && h(7)~=0) || (h(6)~=0 && h(8)~=0)
                        % Surface 3: eta=1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(G1(i,1),1);
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xxi)^2+(Mesh.Elements(j).Shape.Yxi)^2);
                            hInt=[h(5);h(6)]*(1-(1+G1(i,1))/2)+(1+G1(i,1))/2*[h(7);h(8)];
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*hInt*js*G1(i,2);
                                end
                            end
                        end
                    end
                    if (h(1)~=0 && h(7)~=0) || (h(2)~=0 && h(8)~=0)
                        % Surface 4: xi=-1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(-1,G1(i,1));
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xeta)^2+(Mesh.Elements(j).Shape.Yeta)^2);
                            hInt=[h(7);h(8)]*(1-(1+G1(i,1))/2)+(1+G1(i,1))/2*[h(1);h(2)];
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*hInt*js*G1(i,2);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function U=returnInterpolatedU(obj,x)
            U=zeros(2,1);
            for i=1:obj.numberOfNodes
                U=U+obj.Nodes(i).sF.getValue(x)*obj.Nodes(i).u;
            end
        end
    end
end


