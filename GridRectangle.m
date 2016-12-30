function [ NN,NEL,X,Y ] = GridRectangle(xdim,ydim,numElx,numEly)
% Function creates a rectangular grid of dim: xdim,ydim using num elements
% listed
n=numElx*numEly;
dx=(xdim)/numElx; % Node Spacing
dy=(ydim)/numEly; % Node Spacing
nN=(numElx+1)*(numEly+1); % Number of Nodes
[X,Y]=meshgrid([0:dx:xdim],[0:dy:ydim]);

% Define Arrays:
dofx=1:2:nN*2;
dofy=2:2:nN*2;
NN= [(1:nN)',reshape(X',[],1),reshape(Y',[],1),dofx',dofy']; % Global Nodes
NEL=zeros(n,4);

for i=1:n % Loops through the elements
   row=ceil(i/numElx);
   column=mod(i-1,(numElx))+1;
   NEL(i,:)=[(row-1)*(numElx+1)+column,(row-1)*(numElx+1)+column+1,(row)*(numElx+1)+column+1,(row)*(numElx+1)+column];
end

end