%% Create Rectangular Mesh
%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(2,2,50,50);

%Set up Essential Boundary:
b1=[-eps,eps,-eps,eps,[0,0]];
b2=[-eps,eps,1-eps,1+eps,[0,0]];
b3=[-eps,eps,0.5-eps,0.5+eps,[0,0]];
BE=Boundary(NN,b1,b2,b3);

% Set up Natural Boundary:
b2=[2-eps,2+eps,-eps,1+eps,[0,-5000]];
BN=Boundary(NN,b2); BN(BN==-Inf)=0;

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

%% Plot Weight Function(Box):
a=1;
xI=[0,0];
[X,Y]=meshgrid(-1:0.1:1);
Z=zeros(size(X,1));
W=Weight(xI,a,0);
for i=1:size(X,1) % Columns
    for j=1:size(Y,2) % Rows
        Z(i,j)=W.w([X(i,j),Y(i,j)]);
    end
end
surf(X,Y,Z)

%% Plot Weight Function(Box):
a=1;
xI=[0,0];
[X,Y]=meshgrid(-1:0.1:1);
Z=zeros(size(X,1));
W=Weight(xI,a,0);
for i=1:size(X,1) % Columns
    for j=1:size(Y,2) % Rows
        z=W.wx([X(i,j),Y(i,j)]);
        Z(i,j)=z(2);
    end
end
surf(X,Y,Z)



%% Test Point Cloud:
a=0.45;
order=1;
N0=(1:150)';
xNodes=rand(length(N0),1)*2;
yNodes=rand(length(xNodes),1)*2;
a=ones(length(xNodes),1)*a;
padding=zeros(length(xNodes),1);
Nodes=[N0,xNodes,yNodes,a,padding];
PointCloud=Cloud(Nodes,order);
%% Plot Shape Function
Z=zeros(size(X,1));
shape=30;
for i=1:size(X,1) % Columns
    for j=1:size(Y,2) % Rows  
        Z(i,j)=PointCloud.Nodes(shape).sF.getValue([X(i,j),Y(i,j)]);
    end
end
surf(X,Y,Z,'facecolor','blue','EdgeColor','blue')
alpha(0.1)
hold on
for i=1:PointCloud.numberOfNodes
    scatter(PointCloud.Nodes(i).cordinates(1),PointCloud.Nodes(i).cordinates(2),'ko','filled')
    hold on
end

scatter(PointCloud.Nodes(shape).cordinates(1),PointCloud.Nodes(shape).cordinates(2),'r+')


%% Test PU:
xTest=[rand()*2,rand()*2];
pu=0;
for i=1:PointCloud.numberOfNodes
   pu=pu+PointCloud.Nodes(i).sF.getValue([xTest(1),xTest(2)]);
end
assert(abs(pu-1)<eps*10)
