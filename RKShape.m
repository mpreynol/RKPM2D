classdef RKShape < handle
    %Class defines a RKPM shape function 
    % Currently implemeted in 2D!!!
    % Written by Mathew Reynolds on Dec 26, 2016
    
    properties
        a; % Dilation Parameter (Double)
        cordinates=[]; % Cordinates (list, Double)
        order; % Order of polynomial basis used
        Cloud; % Cloud Data Handle
        weight; % Weight Data Handle (from This Node)
        M; % Moment Matrix
        Minv; % Stored value of matrix Invese
        Mdx; % Moment Matrix Derivative.
        Mdy; % Moment Matrix Derivative.
        Minvdx; % Inverse Derivative.
        Minvdy; % Inverse Derivative.
    end
    
    methods
        % Constructor:
        function obj = RKShape(cordinates,a,order,Cloud,weight)
            obj.cordinates=cordinates;
            obj.a=a;
            obj.order=order;
            obj.Cloud=Cloud;
            obj.weight=weight;
        end
                     
        % Basis Vectors
        function H=H(obj,x)
            H=zeros(obj.order*2+1,1);
            H(1)=1;
            counter=2;
            for i=1:obj.order
                H(counter)=(x(1)-obj.cordinates(1))^i;
                H(counter+1)=(x(2)-obj.cordinates(2))^i;
                counter=counter+2;
            end
        end
        function H=Hd(obj,x)
            % order*2+1 rows + SD rows
            % the first column is H,x and the second colum is H,y
            H=zeros(obj.order*2+1,2);
            H(1,:)=0;
            counter=2;
            for i=1:obj.order
                H(counter,1)=i*(x(1)-obj.cordinates(1))^(i-1);
                H(counter+1,1)=0;
                H(counter,2)=0;
                H(counter+1,2)=i*(x(2)-obj.cordinates(2))^(i-1);
                counter=counter+2;
            end
        end
        
        % Processing for Value:
        % Define Moment Matrix:
        function setMoment(obj,x,wAQ)
            Mi=zeros(obj.order*2+1);
            for i=1:(obj.Cloud.numberOfNodes) % Loop through other Nodes in the Mesh
                if wAQ(1,i)~=0
                    H=obj.Cloud.Nodes(i).sF.H(x);
                    Mi=Mi+H*H'*wAQ(1,i);
                end
            end
            obj.M=Mi;
            if sum(isnan(obj.M))>0
                obj.M=zeros(obj.order*2+1);
                obj.Minv=obj.M;
            else
               obj.Minv=inv(obj.M); 
            end
        end
        
        % Processing for Value:
        % Define Moment Matrix:
        function setMomentVanilla(obj,x)
            Mi=zeros(obj.order*2+1);
            for i=1:(obj.Cloud.numberOfNodes) % Loop through other Nodes in the Mesh
                if (obj.Cloud.Nodes(i).cordinates(1)-x(1))<=(obj.Cloud.Nodes(i).a+obj.a) && (obj.Cloud.Nodes(i).cordinates(2)-x(2))<=(obj.Cloud.Nodes(i).a+obj.a) % Only use Nodes from inside dilation of current evaluatioin point
                    H=obj.Cloud.Nodes(i).sF.H(x);
                    Mi=Mi+H*H'*obj.Cloud.Nodes(i).weight.w(x);
                end
            end
            obj.M=Mi;
            if sum(isnan(obj.M))>0
                obj.M=zeros(obj.order*2+1);
                obj.Minv=obj.M;
            else
               obj.Minv=inv(obj.M); 
            end
            
        end 
        
        % Define Value
        function v=getValue(obj,x,wAQ,nodeNo)
            if (prod(x==obj.cordinates) && obj.weight.singular)
                v=1;
            else
                if nargin>2
                    w=wAQ(1,nodeNo);
                    obj.setMoment(x,wAQ);
                else
                    w=obj.weight.w(x);
                    obj.setMomentVanilla(x)
                end
                v=obj.H(obj.cordinates)'*obj.Minv*obj.H(x)*w;
            end      
        end
        
       
        % Processing for Derivative 
        function setMomentdx(obj,x,wAQ)
            Midx=zeros(obj.order*2+1);
            Midy=Midx;
            for i=1:(obj.Cloud.numberOfNodes)
                if wAQ(1,i)~=0
                    % Define Preliminary Values:
                    W=wAQ(1,i);
                    H=obj.Cloud.Nodes(i).sF.H(x);
                    DH=obj.Cloud.Nodes(i).sF.Hd(x);
                    DW=wAQ(2:3,i);
                    % Compute Moment Matrix:
                    Midx=Midx+DH(:,1)*H'*W+H*DH(:,1)'*W+H*H'*DW(1);
                    Midy=Midy+DH(:,2)*H'*W+H*DH(:,2)'*W+H*H'*DW(2);
                end
            end     
            obj.Mdx=Midx;
            obj.Mdy=Midy;
            
            % Compute Inverse Moment Matrix
            obj.Minvdx=-obj.Minv*obj.Mdx*obj.Minv;
            obj.Minvdy=-obj.Minv*obj.Mdy*obj.Minv;
        end
        
        % Define Derivative:
        function vDx=getValueDx(obj,x,wAQ,nodeNo) % Returns [dx;dy]
            vDx=zeros(length(x),1);
            obj.setMoment(x,wAQ);
            obj.setMomentdx(x,wAQ);
            H0=obj.H(obj.cordinates);
            W=wAQ(1,nodeNo);
            H=obj.H(x);
            DH=obj.Hd(x);
            DW=wAQ(2:3,nodeNo);
            vDx(1)=H0'*(obj.Minvdx*H*W+obj.Minv*DH(:,1)*W+obj.Minv*H*DW(1));
            vDx(2)=H0'*(obj.Minvdy*H*W+obj.Minv*DH(:,2)*W+obj.Minv*H*DW(2));
        end
        
    end
    
end

