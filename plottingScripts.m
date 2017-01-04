%% Plot Weight Function(Box):
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





%% Plot Shsape Functions
    
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
