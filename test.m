%% Create Rectangular Mesh
%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(2,2,5,5);



Mesh=QuadMesh(NN,NEL,2);

%% Plot Mesh:
MeshPlot.plotOriginal(Mesh.Elements)

%% Plot Integration Points
for i=1:Mesh.noElements
   for j=1:size(Mesh.Elements(i).G2,1)
       Cords=Mesh.Elements(i).getIntCord(j);
       plot(Cords(1),Cords(2),'b.')
       hold on
   end
end

% %% Plot Weight Function(Box):
% a=1;
% xI=[0,0];
% [X,Y]=meshgrid(-1:0.1:1);
% Z=zeros(size(X,1));
% W=Weight(xI,a,0);
% for i=1:size(X,1) % Columns
%     for j=1:size(Y,2) % Rows
%         Z(i,j)=W.w([X(i,j),Y(i,j)]);
%     end
% end
% surf(X,Y,Z)
% 
% %% Plot Weight Function Derivative(Box):
% a=1;
% xI=[0,0];
% [X,Y]=meshgrid(-1:0.1:1);
% Z=zeros(size(X,1));
% W=Weight(xI,a,0);
% for i=1:size(X,1) % Columns
%     for j=1:size(Y,2) % Rows
%         z=W.wx([X(i,j),Y(i,j)]);
%         Z(i,j)=z(2);
%     end
% end
% surf(X,Y,Z)

%% Generate Point Cloud:
[NNp] = GridRectangle(2,2,5,5);
%Set up Essential Boundary:
b1=[-eps,eps,-eps,eps,[0,0]];
b2=[-eps,eps,1-eps,1+eps,[0,0]];
b3=[-eps,eps,0.5-eps,0.5+eps,[0,0]];
BE=Boundary(NNp,b1,b2,b3);

% Set up Natural Boundary:
b2=[2-eps,2+eps,-eps,1+eps,[0,-5000]];
BN=Boundary(NNp,b2); BN(BN==-Inf)=0;
%Points=Mesh.createPoints(20);

%% Test Point Cloud:
a=0.5;
order=1;
xNodes=NNp(:,2);
NO=(1:length(xNodes))';
yNodes=NNp(:,3);
A=ones(size(NN,1),1)*a;
padding=zeros(length(A),1);
Nodes=[NO,xNodes,yNodes,ones(length(NO),1)*a,zeros(length(NO),1)];
PointCloud=Cloud(Nodes,order);

%% Plot Mesh
    for i=1:PointCloud.numberOfNodes
        refresh
        scatter(PointCloud.nodalData(i,2),PointCloud.nodalData(i,3),'ko','filled')
        hold on
    end
% %% Plot Shape Function
% Z=zeros(size(X,1));
% shape=20;
% for i=1:size(X,1) % Columns
%     for j=1:size(Y,2) % Rows  
%         Z(i,j)=PointCloud.Nodes(shape).sF.getValue([X(i,j),Y(i,j)]);
%     end
% end
% surf(X,Y,Z,'facecolor','blue','EdgeColor','blue')
% alpha(0.1)
% hold on
% for i=1:PointCloud.numberOfNodes
%     scatter(PointCloud.Nodes(i).cordinates(1),PointCloud.Nodes(i).cordinates(2),'ko','filled')
%     hold on
% end
% 
% scatter(PointCloud.Nodes(shape).cordinates(1),PointCloud.Nodes(shape).cordinates(2),'r+')
% 
% %% Test PU:
% xTest=[rand()*2,rand()*2];
% pu=0;
% for i=1:PointCloud.numberOfNodes
%    pu=pu+PointCloud.Nodes(i).sF.getValue([xTest(1),xTest(2)]);
% end
% assert(abs(pu-1)<eps*10)
% 
% %% Plot Shape Function Derivative
% Z1=zeros(size(X,1));
% Z2=zeros(size(X,1));
% for i=1:size(X,1) % Columns
%     for j=1:size(Y,2) % Rows  
%         z=PointCloud.Nodes(shape).sF.getValueDx([X(i,j),Y(i,j)]);
%         Z1(i,j)=z(1); Z2(i,j)=z(2);
%     end
% end
% % Plots!
% subplot(1,2,1)
% surf(X,Y,Z1,'facecolor','blue','EdgeColor','blue')
% alpha(0.1)
% hold on
% for i=1:PointCloud.numberOfNodes
%     scatter(PointCloud.Nodes(i).cordinates(1),PointCloud.Nodes(i).cordinates(2),'ko','filled')
%     hold on
% end
% scatter(PointCloud.Nodes(shape).cordinates(1),PointCloud.Nodes(shape).cordinates(2),'r+')
% 
% subplot(1,2,2)
% surf(X,Y,Z2,'facecolor','blue','EdgeColor','blue')
% alpha(0.1)
% hold on
% for i=1:PointCloud.numberOfNodes
%     scatter(PointCloud.Nodes(i).cordinates(1),PointCloud.Nodes(i).cordinates(2),'ko','filled')
%     hold on
% end
% scatter(PointCloud.Nodes(shape).cordinates(1),PointCloud.Nodes(shape).cordinates(2),'r+')

%% Integration Loop:
C=Constit(1E6,0.2,'Plane Stress').C;
n=PointCloud.numberOfNodes;
K=zeros(n*2);
F=zeros(n*2,1);
for k=1:Mesh.noElements % Loop through Elements
    k
    for l=1:size(Mesh.Elements(k).G2,1) % Loop through Quadrature Points
        Cords=Mesh.Elements(k).getIntCord(l); x=Cords(1); y=Cords(2); Jw=Cords(3); % Jw = Jacobian * Quadrature Weight
        for a=1:n % Loop over Shape Functions
            for b=1:n % Loop over Shape Functions
                if b>=a
                    CordsA=PointCloud.Nodes(a).cordinates; CordsB=PointCloud.Nodes(b).cordinates;
                    if ((CordsA(1)-CordsB(1))^2+(CordsA(2)-CordsB(2))^2)<PointCloud.Nodes(a).a^2
                        NDa=PointCloud.Nodes(a).sF.getValueDx([x;y]);
                        NDb=PointCloud.Nodes(b).sF.getValueDx([x;y]);
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
            F(2*a-1)=0;%F(2*a-1)+PointCloud.Nodes(a).sF.getValue([x;y]);
            F(2*a)=F(2*a)+PointCloud.Nodes(a).sF.getValue([x;y])*-1;
        end
    end
end
K=triu(K)+tril(K',-1);
%% Solve System:
uUnknown=ones(n,1)*-inf; uUnknown(1)=0; logical=(uUnknown==-inf);
Kr=K(logical,logical); Fr=F(logical);
dr=Kr\Fr;
d=[0;dr];

%% Interpolate the Solution
Result=zeros(length(xI),1);
for g=1:length(xI)
    x=xI(g);
    for i=1:n
       Result(g)=Result(g)+Shapes(i).getValue(x)*d(i);
    end
end
plot(xI,Result)
