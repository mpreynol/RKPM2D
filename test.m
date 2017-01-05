%% Create Rectangular Mesh for Quadrature
%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(2,2,20,20);
Mesh=QuadMesh(NN,NEL,2);

%% Plot Mesh:
MeshPlot.plotOriginal(Mesh.Elements)
MeshPlot.plotQuadraturePoints(Mesh)

%% Generate Point Cloud:
[NNp] = GridRectangle(2,2,6,6);
%Set up Essential Boundary:
b1=[-eps,eps,-eps,2+eps,[0,0]];
b2=[-eps,eps,1-eps,1+eps,[0,0]];
b3=[-eps,eps,2-eps,2+eps,[0,0]];
BE=Boundary(NNp,b1);

% Set up Natural Boundary:
b2=[2-eps,2+eps,-eps,1+eps,[0,-5000]];
BN=Boundary(NNp,b2); BN(BN==-Inf)=0;
%Points=Mesh.createPoints(20);

%% Define Point Cloud Object:
a=sqrt(2*(0.33333334)^2)*1.01;
order=1;
xNodes=NNp(:,2);
yNodes=NNp(:,3);
Nodes=[(1:length(xNodes))',xNodes,yNodes,ones(length(xNodes),1)*a,reshape(BE,2,[])'];
PointCloud=Cloud(Nodes,order);
PointCloud.plotCloud()

%% Build Arrays:
C=Constit(1,0,'Plane Stress').C;
Q=[0;-.1];
[K,F]=PointCloud.integrateDomain(C,Q,Mesh);


%% Solve System:
L=BE==-inf; % Indexes of unknown equations
Kr=K(L,L); Br=BE(~L); fr=F(L); KRHS=K(L,~L); RHS=fr-KRHS*Br;
ur=Kr\RHS;
u=PointCloud.reAssembleUnknowns(ur,BE);

% Populate solution back into PointCloud Collection:
PointCloud.parseSolution(u);

%% Plot The Solution:
PointCloud.plotCloud();
PointCloud.plotCloudU(10)

%% Test Points
xs=0:0.1:2;
for q=1:length(xs)
    ui=PointCloud.returnInterpolatedU([xs(q);0]);
    plot(xs(q),ui(2),'.');
    hold on
end

%%
PointCloud.returnInterpolatedU([2;0])
