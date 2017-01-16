%% Create Rectangular Mesh for Quadrature
%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(20,1,40,4);
orderInt=3;
threads=4;
%Mesh=QuadMesh(NN,NEL,orderInt);
SplitMesh=Domains(NN,NEL,orderInt,threads);

% Set up Natural Boundary:
b2=[20-eps,20+eps,-eps,1+eps,[0,1]];
BN=Boundary(NN,b2); BN(BN==-Inf)=0;

%% Plot Mesh:
%MeshPlot.plotOriginalDomains(SplitMesh);
%MeshPlot.plotQuadraturePointsDomains(SplitMesh);

%% Generate Uniform Point Cloud:
[NNp] = GridRectangle(20,1,40,4);
%Set up Essential Boundary:
b1=[-eps,eps,0.5-eps,0.5+eps,[0,0]];
b2=[-eps,eps,0-eps,0+eps,[0,-Inf]];
b3=[-eps,eps,1-eps,1+eps,[0,-Inf]];
BE=Boundary(NNp,b1,b2,b3);
xNodes=NNp(:,2);
yNodes=NNp(:,3);
Nodes=[(1:length(xNodes))',xNodes,yNodes,reshape(BE,2,[])'];
%PointCloud=Cloud(Nodes,1,1);
PointCloud2=Cloud2(Nodes,1,1);
%PointCloud.plotCloud()

% %% Define Random Point Cloud Object:
% Points=Mesh.createPoints(60,[0,0;0,0.5;0,1;1,1;1,0]);
% %Set up Essential Boundary:
% b1=[-eps,eps,0.5-eps,0.5+eps,[0,0]];
% b2=[-eps,eps,0-eps,0+eps,[0,0]];
% b3=[-eps,eps,1-eps,1+eps,[0,0]];
% BE=Boundary([(1:size(Points,1))',Points,reshape((1:size(Points,1)*2),2,[])'],b1,b2,b3);
% 
% xNodes=Points(:,1);
% yNodes=Points(:,2);
% Nodes=[(1:length(xNodes))',xNodes,yNodes,reshape(BE,2,[])'];
% PointCloud=Cloud(Nodes,1,2);
% PointCloud.plotCloud()

%% Perform Checks

%PointCloud.checkPU(Mesh)
%PointCloud.checkNU(Mesh)
%PointCloud.checkLinearConsistency(Mesh)
%% Precalculate
SplitMesh.weightatQuad(PointCloud2);
SplitMesh.RKPMatQuad(PointCloud2); % Store Values of Shape Functions W/ decomposition
%Mesh.WeightatQuadrature(PointCloud);
%[RKv,RKdx,DRKdy]=Mesh.RKPMatQuadrature(PointCloud); 
C=Constit(1E6,0.3,'Plane Stress').C;
%% Integration
Q=[0;0];
[K,Fb]=PointCloud2.integrateDomain(C,Q,SplitMesh); % Build Arrays
%[K,Fb]=PointCloud.integrateDomain(C,Q,Mesh,RKv,RKdx,DRKdy); % Build Arrays
%Fh=PointCloud.integrateExactBoundary(Mesh,BN,@parabolicStress); sum(Fh);
Fh=PointCloud2.integrateExactBoundary(SplitMesh,BN,@parabolicStress); sum(Fh);
F=Fb+Fh;
%% Solve System:
L=BE==-inf; % Indexes of unknown equations
Kr=K(L,L); Br=BE(~L); fr=F(L); KRHS=K(L,~L); RHS=fr-KRHS*Br;
ur=Kr\RHS;
u=PointCloud2.reAssembleUnknowns(ur,BE);
% Populate solution back into PointCloud Collection:
PointCloud2.parseSolution(u);

%% Plot The Solution:
%PointCloud.plotCloud();
%d=PointCloud.plotCloudU(1);

% %% Test Points
% xs=0:0.1:1;
% for q=1:length(xs)
%     ui=PointCloud.returnInterpolatedU([xs(q);1]);
%     plot(xs(q),ui(2),'k.');
%     hold on
% end
% 

%%
PointCloud2.returnInterpolatedU([20;0])
