%% Create Rectangular Mesh for Quadrature
%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(2,2,60,60);
Mesh=QuadMesh(NN,NEL,2);

% Set up Natural Boundary:
b2=[2-eps,2+eps,-eps,1+eps,[-0.1,0]];
BN=Boundary(NN,b2); BN(BN==-Inf)=0;

%% Plot Mesh:
MeshPlot.plotOriginal(Mesh.Elements)
MeshPlot.plotQuadraturePoints(Mesh)

%% Generate Point Cloud:
[NNp] = GridRectangle(2,2,6,6);
%Set up Essential Boundary:
b1=[-eps,eps,1-eps,1+eps,[0,0]];
b2=[-eps,eps,-eps,2+eps,[0,0]];
b3=[-eps,eps,1-eps,1+eps,[0,0]];
BE=Boundary(NNp,b2,b1);

%% Define Point Cloud Object:
a=sqrt(0.33334^2+0.333334^2)*1.01*2;
%Points=Mesh.createPoints(20);
order=1;
xNodes=NNp(:,2);
yNodes=NNp(:,3);
Nodes=[(1:length(xNodes))',xNodes,yNodes,ones(length(xNodes),1)*a,reshape(BE,2,[])'];
PointCloud=Cloud(Nodes,order);
PointCloud.plotCloud()

%% Build Arrays:
C=Constit(1,0,'Plane Strain').C;
Mesh.RKPMatQuadrature(PointCloud); % Store Values of Shape Functions 
Q=[0;0];
[K,Fb]=PointCloud.integrateDomain(C,Q,Mesh); % Build Arrays
Fh=PointCloud.integrateBoundary(Mesh,BN); sum(Fh)
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
PointCloud.plotCloudU(1)

%% Test Points
xs=0:0.1:2;
for q=1:length(xs)
    ui=PointCloud.returnInterpolatedU([xs(q);1]);
    plot(xs(q),ui(1),'k.');
    hold on
end

%%
PointCloud.returnInterpolatedU([2;0])
