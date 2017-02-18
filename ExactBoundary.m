function [B] = ExactBoundary(NN,varargin)
% Boundary Function takes in Nodal Array and a collection of boundary
% values builds an array NN elements long specifying the assigned
% boundary value or -inf if DOF is free.
% Boundary is set up for 2D right now

% Assumptions
%a. Boundary is the closed set
%b. NN is formatted as per "GridSquare" function
%c. To specify an edge specify a min and max correspoding a thin threshold
D=2;
I=D^3/12;
L=10;
E=21.1E6;
P=-20000;
mu=0.3;


% Syntax for Boundary Specification:
% [xmin, xmax, ymin,ymax,[value]]
B=ones(size(NN,1)*2,1)*-inf; % Length of Boundary Vector is the no DOF
for i=1:nargin-1
    infoArray=varargin{i};
    for j=1:size(NN,1)
        if (NN(j,2)>=infoArray(1) && NN(j,2)<=infoArray(2) && NN(j,3)>=infoArray(3) && NN(j,3)<=infoArray(4))
            B(NN(j,4))=-P.*(NN(j,3)-1)/(6*E*I).*((2+mu).*((NN(j,3)-1).^2-D^2/4));
            B(NN(j,5))=P/(6*E*I)*(3*mu*(NN(j,3)-1).^2.*(L));
        end
    end
end
end