classdef Weight < handle
    %WEIGHT is a class of statics methods for the weighting function and
    %Derivative
    
    properties
        a; % Dilation Parameter
        s; % Weight Function Home location
        p=2; % Order of the Singularity
        singular=false; % Boolean Flag whether or not weight function is singulae
        F; % Value of the singularity Function
        DF; % Value of the singularity Function Derivative
    end
    
    methods
        % Constructor
        function obj=Weight(s,a,varargin)
            obj.s=s;
            obj.a=a;
            if nargin>2
                obj.singular=varargin{1};
            end
        end
        
        % Singular Function 
        function f=f(obj,x)
           f=(((x(1)-obj.s(1))/obj.a)^2+((x(2)-obj.s(2))/obj.a)^2)^obj.p;
        end
        function df=df(obj,x)
            df(1)=2*obj.p*(((x(1)-obj.s(1))/obj.a))^(2*obj.p-1)*(1/obj.a);
            df(2)=2*obj.p*(((x(2)-obj.s(2))/obj.a))^(2*obj.p-1)*(1/obj.a);
        end
        
        function w = w(obj,x) % x,s is an array of length SD
            w=prod(obj.weights(x));
        end
        
        function wx = wx(obj,x) % x,s is an array of length SD
            weights=obj.weights(x);
            derivative=obj.derivative(x);
            wx=zeros(length(x),1);
            temp=1;
            for i=1:length(x)
                for j=1:length(x)
                    if j~=i
                        temp=temp*weights(j);
                    end
                end
                wx(i)=derivative(i)*temp;
            end
        end
        
        % Weight Function
        function weights = weights(obj,x)   
            z=zeros(length(x),1);
            weights=zeros(length(x),1);
            for i=1:length(x);
                z(i)=abs(obj.s(i)-x(i))/obj.a;   
                if z(i)>=1
                    weights(i)=0;
                elseif z(i)<1/2
                    weights(i)=2/3-4*z(i)^2+4*z(i)^3;
                else
                    weights(i)=4/3-4*z(i)+4*z(i)^2-4/3*z(i)^3;
                end
                if obj.singular
                    obj.F=obj.f(x);
                    weights(i)=weights(i)/obj.F;
                end
            end
        end
        
        % Weight Function Derivative
        function wx=derivative(obj,x)
            z=zeros(length(x),1);
            dz=zeros(length(x),1);
            wx=zeros(length(x),1);  
            for i=1:length(x)  
                z(i)=abs(obj.s(i)-x(i))/obj.a;
                dz(i)=(obj.s(i)-x(i))/abs(obj.s(i)-x(i))/obj.a;
                if z(i)<eps && z(i)>-eps
                    dz(i)=0;
                end
                if z(i)>=1
                    wx(i)=0;
                elseif z(i)<1/2
                    wx(i)=-8*z(i)*dz(i)+12*z(i).^2*dz(i);
                else
                    wx(i)=-4*dz(i)+8*z(i)*dz(i)-4*z(i).^2*dz(i);
                end
                if obj.singular
                    obj.F=obj.f(x);
                    obj.DF=obj.df(x);
                    wx(i)=(obj.F*wx(i)-obj.w(x)*obj.DF(i))/(obj.F^2);
                end
            end
        end
    end
    
end