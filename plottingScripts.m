% Plot Weight Function(Box):
a=1;
xI=[0,0];
[X,Y]=meshgrid(-1:0.1:1);
Z=zeros(size(X,1));
W=Weight(xI,a,0);
for i=1:size(X,1) % Columns
    for j=1:size(Y,2) % Rows
        Z(i,j)=W.w([X(i,j);Y(i,j)]);
    end
end
surf(X,Y,Z)

%% Point on Shape Function Derivative:
a=1;
xI=[0,0];
W=Weight(xI,a,0);
z=W.wx([-0.5;-0.5])
    
%% Plot Weight Function Derivative(Box):
a=1;
xI=[0,0];
[X,Y]=meshgrid(-1:0.1:1);
Z=zeros(size(X,1));
W=Weight(xI,a,0);
for i=1:size(X,1) % Columns
    for j=1:size(Y,2) % Rows
        z=W.wx([X(i,j);Y(i,j)]);
        Z(i,j)=z(2);
    end
end
surf(X,Y,Z)


%% Plot Shape Function
figure(1)
shape=11;
Z=zeros(size(X,1));
for i=1:size(X,1) % Columns
    for j=1:size(Y,2) % Rows  
        Z(i,j)=PointCloud.Nodes(shape).sF.getValue([X(i,j);Y(i,j)]);
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
xTest=Mesh.Elements(1).getIntCord(1);
pu=0;
for i=1:PointCloud.numberOfNodes
   pu=pu+PointCloud.Nodes(i).sF.getValue([xTest(1);xTest(2)]);
end
pu-1

%% Test NU:

xTest=Mesh.Elements(1).getIntCord(1);
nu1=0;
nu2=0;
for i=1:PointCloud.numberOfNodes
    temp=PointCloud.Nodes(i).sF.getValueDx([xTest(1);xTest(2)]);
    nu1=nu1+temp(1);
    nu2=nu2+temp(2);
end
nu1
nu2


%% Point on Shape Function Deriviative
shape=1;
T=Mesh.Elements(1).getIntCord(1);
x=[T(1);T(2)];
value=PointCloud.Nodes(shape).sF.getValue(x)
derivative=PointCloud.Nodes(shape).sF.getValueDx(x)

%% Plot Shape Function Derivative
figure(2)
Z1=zeros(size(X,1));
Z2=zeros(size(X,1));
for i=1:size(X,1) % Columns
    for j=1:size(Y,2) % Rows  
        if 1
            z=PointCloud.Nodes(shape).sF.getValueDx([X(i,j),Y(i,j)]);
            Z1(i,j)=z(1); Z2(i,j)=z(2);
        else
            % Pass
        end
    end
end
% Plots!
subplot(1,2,1)
surf(X,Y,Z1,'facecolor','blue','EdgeColor','blue')
alpha(0.1)
hold on
for i=1:PointCloud.numberOfNodes
    scatter(PointCloud.Nodes(i).cordinates(1),PointCloud.Nodes(i).cordinates(2),'ko','filled')
    hold on
end
scatter(PointCloud.Nodes(shape).cordinates(1),PointCloud.Nodes(shape).cordinates(2),'r+')

subplot(1,2,2)
surf(X,Y,Z2,'facecolor','blue','EdgeColor','blue')
alpha(0.1)
hold on
for i=1:PointCloud.numberOfNodes
    scatter(PointCloud.Nodes(i).cordinates(1),PointCloud.Nodes(i).cordinates(2),'ko','filled')
    hold on
end
scatter(PointCloud.Nodes(shape).cordinates(1),PointCloud.Nodes(shape).cordinates(2),'r+')
