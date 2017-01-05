Nodes=[0,0;1,0;2,0;0,1;1,1;2,1;0,2;1,2;2,2];
x=[1;0.5];
a=1.5;
%% Compute Weight Function Values
for i=1:9
    W=Weight([Nodes(i,:)'],a);
    W.weights(x);
    W.w(x);
    %Node=RKNode(i,[Nodes(i,:)'],a,1,[-Inf,-Inf],Cloud)
end
%% Create Cloud Object
PointCloud=Cloud([(1:9)',Nodes(:,1),Nodes(:,2),ones(9,1)*a,ones(9,1)*-Inf,ones(9,1)*-Inf],1);

%% Test PU:
pu=0;
for i=1:PointCloud.numberOfNodes
    pu=pu+PointCloud.Nodes(i).sF.getValue(x);
end
