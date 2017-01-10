%% Create Rectangular Mesh for Quadrature
%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(1,1,6,6);
Mesh=QuadMesh(NN,NEL,2);

% Set up Natural Boundary:
b2=[1-eps,1+eps,-eps,1+eps,[0.5,-0.1]];
BN=Boundary(NN,b2); BN(BN==-Inf)=0;

%% Plot Mesh:
MeshPlot.plotOriginal(Mesh.Elements)
MeshPlot.plotQuadraturePoints(Mesh)

%% Generate Uniform Point Cloud:
[NNp] = GridRectangle(1,1,6,6);
%Set up Essential Boundary:
b1=[-eps,eps,0.5-eps,0.5+eps,[0,0]];
b2=[-eps,eps,0-eps,0+eps,[0,0]];
b3=[-eps,eps,1-eps,1+eps,[0,0]];
BE=Boundary(NNp,b1,b2,b3);
xNodes=NNp(:,2);
yNodes=NNp(:,3);
Nodes=[(1:length(xNodes))',xNodes,yNodes,reshape(BE,2,[])'];
PointCloud=Cloud(Nodes,1,2);
PointCloud.plotCloud()

%% Define Random Point Cloud Object:
Points=Mesh.createPoints(60,[0,0;0,0.5;0,1;1,1;1,0]);
%Set up Essential Boundary:
b1=[-eps,eps,0.5-eps,0.5+eps,[0,0]];
b2=[-eps,eps,0-eps,0+eps,[0,0]];
b3=[-eps,eps,1-eps,1+eps,[0,0]];
BE=Boundary([(1:size(Points,1))',Points,reshape((1:size(Points,1)*2),2,[])'],b1,b2,b3);

xNodes=Points(:,1);
yNodes=Points(:,2);
Nodes=[(1:length(xNodes))',xNodes,yNodes,reshape(BE,2,[])'];
PointCloud=Cloud(Nodes,1,2);
PointCloud.plotCloud()

%% Perform Checks
PointCloud.checkPU(Mesh)
PointCloud.checkNU(Mesh)
PointCloud.checkLinearConsistency(Mesh)

%% Build Arrays:
C=Constit(1,0.2,'Plane Strain').C;
Mesh.RKPMatQuadrature(PointCloud); % Store Values of Shape Functions 
Q=[0;0];
[K,Fb]=PointCloud.integrateDomain(C,Q,Mesh); % Build Arrays
Fh=PointCloud.integrateBoundary(Mesh,BN); sum(Fh);
F=Fb+Fh;
%% Solve System:
L=BE==-inf; % Indexes of unknown equations
Kr=K(L,L); Br=BE(~L); fr=F(L); KRHS=K(L,~L); RHS=fr-KRHS*Br;
ur=Kr\RHS;
u=PointCloud.reAssembleUnknowns(ur,BE);
% Populate solution back into PointCloud Collection:
PointCloud.parseSolution(u);

%% Plot The Solution:
PointCloud.plotCloud();
d=PointCloud.plotCloudU(0.2);

%% Test Points
xs=0:0.1:1;
for q=1:length(xs)
    ui=PointCloud.returnInterpolatedU([xs(q);1]);
    plot(xs(q),ui(2),'k.');
    hold on
end

%%
 PointCloud.returnInterpolatedU([1;0])
