classdef Domains < handle
    
    
    properties
        noDomains; % Number of Domains
        domainCollection; % Collection for Domains
        orderInt; % level of quadrature for mesh
        NN; % Array of Grid Points
        NEL; % Element Connectivity
        noElements; % No. for elements
        commonDomainSize; %Size of common Domain
        unCommonDomainSize; % size of uncommon Domain
    end
    
    methods
        function obj = Domains(NN,NEL,orderInt,noDomains)
           obj.domainCollection=QuadMesh.empty(noDomains,0); % Declare Empty Collection 
           obj.noDomains=noDomains;
           obj.orderInt=orderInt;
           obj.NN=NN;
           obj.NEL=NEL;
           obj.noElements=size(NEL,1);
           obj.commonDomainSize=floor(obj.noElements/noDomains);
           obj.unCommonDomainSize=rem(obj.noElements,noDomains)+obj.commonDomainSize;
           obj.defineDomains();
        end
        
        function defineDomains(obj)
           a=1;
           for i=1:obj.noDomains
               if i==obj.noDomains
                   b=size(obj.NEL,1);
               else
                b=a+obj.commonDomainSize-1;
               end
               localNEL=obj.NEL(a:b,:);
               obj.domainCollection(i)=QuadMesh(obj.NN,localNEL,obj.orderInt);
               a=b+1;
           end
        end
        
        function RKPMatQuad(obj,Cloud)
            for i=1:obj.noDomains
                obj.domainCollection(i).RKPMatQuadrature(Cloud);
            end
        end
        
        function weightatQuad(obj,Cloud)
            for i=1:obj.noDomains
                obj.domainCollection(i).WeightatQuadrature(Cloud);
            end
        end
    end
    
end

