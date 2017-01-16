classdef Weight < handle
    %WEIGHT is a class of statics methods for the weighting function and
    %Derivative
    properties
        a; % Dilation Parameter
        s; % Weight Function Home location
        p=1/2; % Order of the Singularity
        singular=false; % Boolean Flag whether or not weight function is singulae
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
       
        % Singular Function (Boxed Domain?)
        function f=f(obj,x)
           f=(((x(1)-obj.s(1))/obj.a)^2+((x(2)-obj.s(2))/obj.a)^2)^obj.p;
        end
        function df=df(obj,x)
            df(1)=obj.p*(((x(1)-obj.s(1))/obj.a)^2+((x(2)-obj.s(2))/obj.a)^2)^(obj.p-1)*2*(x(1)-obj.s(1))/(obj.a^2);
            df(2)=obj.p*(((x(1)-obj.s(1))/obj.a)^2+((x(2)-obj.s(2))/obj.a)^2)^(obj.p-1)*2*(x(2)-obj.s(2))/(obj.a^2);
        end        
        function w = w(obj,x,weights) % x,s is an array of length SD
            if nargin>2
                w=weights;
            else
                w=obj.weights(x);
            end
            w=w(1)*w(2);
            if obj.singular
                    F=obj.f(x);
                    w=w/F;
            end
        end        
        function wx = wx(obj,x,weights) % x,s is an array of length SD
            if nargin>2
                %weights;
            else
                weights=obj.weights(x);
            end
            derivative=obj.derivative(x);
            wx=zeros(length(x),1);
            F=obj.f(x);
            DF=obj.df(x);
            wx(1)=derivative(1)*weights(2);
            wx(2)=weights(1)*derivative(2);
            if obj.singular
                wx(1)=((F*derivative(1)-weights(1)*DF(1))/F^2)*weights(2);
                wx(2)=((F*derivative(2)-weights(2)*DF(2))/F^2)*weights(1);
            end
        end       
        % Weight Function
        function weights = weights(obj,x)   
            z=[0;0];
            weights=[0;0];
            for i=1:length(x);
                z(i)=abs(x(i)-obj.s(i))/obj.a;   
                if z(i)>=1
                    weights(i)=0;
                elseif z(i)<1/2
                    weights(i)=2/3-4*z(i)^2+4*z(i)^3;
                else
                    weights(i)=4/3-4*z(i)+4*z(i)^2-4/3*z(i)^3;
                end
            end
        end        
        % Weight Function Derivative
        function wx=derivative(obj,x)
            z=[0;0];
            dz=[0;0];
            wx=[0;0];  
            for i=1:length(x)  
                z(i)=abs(x(i)-obj.s(i))/obj.a;   
                dz(i)=(x(i)-obj.s(i))/abs(x(i)-obj.s(i))/obj.a;
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
            end
        end
    end   
end