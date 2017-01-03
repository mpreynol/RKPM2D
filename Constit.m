classdef Constit<handle
    %CONSTIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lambda; % Lames Constant
        mu; % Lames Constant
        nu; % Poissons Ratio
        E; % Youngs Modulus
        alpha; % Constiutive Switching property
        token; % Plane Strain or Plane Stress Not case or space sensitive
        C; % Constituitive Tensor
    end
    
    methods
        % Overloaded Constructor
        function obj = Constit(E,nu,token)
            obj.E=E;
            obj.nu=nu;
            obj.token=lower(strtrim(token));
            obj.lambda=nu*E/((1+nu)*(1-2*nu));
            obj.mu=E/(2*(1+nu));
            obj.alpha=obj.lambda/(obj.lambda+2*obj.mu);
            obj.C=zeros(3,3);
            obj.setTensor()
        end
        
        %Method that returns 3x3 constiutive tensor
        function setTensor(obj)
            obj.C(1,:)=[obj.lambda + 2*obj.mu, obj.lambda, 0];
            obj.C(2,:)=[obj.lambda, obj.lambda + 2*obj.mu, 0];
            obj.C(3,:) = [0,0,obj.mu];
            switch obj.token
                case 'planestrain'
                    obj.C=obj.C; % Unmodified Strain Tensor
                case 'planestress'
                    shift=[-obj.alpha*obj.mu,-obj.alpha*obj.mu,0;-obj.alpha*obj.mu,-obj.alpha*obj.mu,0;0,0,0];
                    obj.C=obj.C + shift;
            end
        end
        
    end
    
end

