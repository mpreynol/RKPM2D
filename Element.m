classdef Element < handle
    %ELEMENT is a general element for quadrature in RKPM
    % Linked with FEM Shape Functions for the quadrature utilities
    
    properties
        nodes=[];
        x=[]; % List of x cordinates for element
        y=[]; % List of y cordinates for element
        orderInt; % Order of integration to develop components
        dof; % DOF List
        G2=[]; % Guassian Integration Array in 2D
        G1=[]; % Guassian Integration Array in 1D
        Shape; % object containing shape functions
        StoredRKPM=[]; % Values of RKPM Shape Functions stored at each integration point size: (1+SD)x(NP)xGP
        % Order of RKPM array is same as G2 Array!!!
        % The third Dimension of the RKPM array handles gauss points
        StoredWeight=[]; % Stored Value of Weight Functions
    end
    
    methods
        % Overloaded constructor with nodal inputs and integration order
        function obj = Element(x,y,nodes,dof,orderInt)
            obj.x=x;
            obj.y=y;
            obj.nodes=nodes;
            obj.dof=dof;
            obj.orderInt=orderInt;
            obj.G2=Quadrature.twoDim(orderInt);
            obj.G1=Quadrature.oneDim(orderInt);
            obj.Shape=quadLinear(x,y); % Initialize Shape Function
        end
        
        function getIntCord=getIntCord(obj,number)
            % Returns in 2D a row vector of: x,y cordinates, integration_weight*jacobian
           parametricCords=[obj.G2(number,1),obj.G2(number,2)]; 
           % Set Shape Functions to get our Mapping:
           obj.Shape.setAll(parametricCords(1),parametricCords(2));
           getIntCord=[obj.Shape.X, obj.Shape.Y, obj.Shape.j*obj.G2(number,3)];
        end
        
        function setStoredRKPM(obj,Value)
            obj.StoredRKPM=Value;
        end
        
        function setStoredWeight(obj,Value)
            %Dimensions: (1+SD)xNPxGP, 3rd dimension depends on gauss order
            obj.StoredWeight=Value;
        end
        
    end
    
end

