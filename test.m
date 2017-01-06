%% Create Rectangular Mesh for Quadrature
%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(2,2,8,8);
Mesh=QuadMesh(NN,NEL,3);

%% Plot Mesh:
MeshPlot.plotOriginal(Mesh.Elements)
MeshPlot.plotQuadraturePoints(Mesh)

%% Generate Point Cloud:
[NNp] = GridRectangle(2,2,8,8);
%Set up Essential Boundary:
b1=[-eps,eps,1-eps,1+eps,[0,0]];
b2=[-eps,eps,-eps,2+eps,[0,-Inf]];
b3=[-eps,eps,2-eps,2+eps,[0,0]];
BE=Boundary(NNp,b1,b2);

% Set up Natural Boundary:
b2=[2-eps,2+eps,-eps,1+eps,[0,-5000]];
BN=Boundary(NNp,b2); BN(BN==-Inf)=0;


%% Define Point Cloud Object:
a=sqrt(0.25^2+0.25^2)*1.01*1;
%Points=Mesh.createPoints(20);
order=1;
xNodes=NNp(:,2);
yNodes=NNp(:,3);
Nodes=[(1:length(xNodes))',xNodes,yNodes,ones(length(xNodes),1)*a,reshape(BE,2,[])'];
PointCloud=Cloud(Nodes,order);
PointCloud.plotCloud()

%% Build Arrays:
C=Constit(1,0.2,'Plane Strain').C;
Mesh.RKPMatQuadrature(PointCloud); % Store Values of Shape Functions 
Q=[-1;0];
[K,F]=PointCloud.integrateDomain(C,Q,Mesh); % Build Arrays


%% Solve System:
L=BE==-inf; % Indexes of unknown equations
Kr=K(L,L); Br=BE(~L); fr=F(L); KRHS=K(L,~L); RHS=fr-KRHS*Br;
ur=Kr\RHS;
u=PointCloud.reAssembleUnknowns(ur,BE);

% Populate solution back into PointCloud Collection:
PointCloud.parseSolution(u);

%% Plot The Solution:
PointCloud.plotCloud();
PointCloud.plotCloudU(0.5)

%% Test Points
xs=0:0.1:2;
for q=1:length(xs)
    ui=PointCloud.returnInterpolatedU([xs(q);0]);
    plot(xs(q),ui(1),'k.');
    hold on
end

%%
PointCloud.returnInterpolatedU([2;0])
