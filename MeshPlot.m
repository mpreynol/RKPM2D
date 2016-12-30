classdef MeshPlot<handle
    %MESHPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function plotOriginal(Mesh)
            for i=1:size(Mesh,2)
                x=Mesh(i).x; x=[x;x(1)];
                y=Mesh(i).y; y=[y;y(1)];
                lh.Color=[0,0,0,0.5];
                plot(x,y,'-',lh)
                hold on
            end
        end
        function plotDeformed(Mesh,scale)
             for i=1:size(Mesh,2)
                x=Mesh(i).x; x=[x;x(1)];
                dofx=Mesh(i).dof(1:2:8);
                dofy=Mesh(i).dof(2:2:8);
                ux=Mesh(i).u(1:2:8); ux=[ux;ux(1)];
                y=Mesh(i).y; y=[y;y(1)];
                uy=Mesh(i).u(2:2:8); uy=[uy;uy(1)];
                plot((x+scale*ux),(y+scale*uy),'-b')
                hold on
            end           
        end
        function plotexx(Mesh)
            maxValue=0;
            minValue=0;
            for i=1:size(Mesh,2)
               if max(Mesh(i).E(1,:))> maxValue
                   maxValue=max(Mesh(i).E(1,:))+eps;
               end
               if min(Mesh(i).E(1,:))< minValue
                   minValue=min(Mesh(i).E(1,:))-eps;
               end               
            end
            for i=1:size(Mesh,2)
                x=Mesh(i).x;
                y=Mesh(i).y;       
                for j=1:4
                    c=Mesh(i).E(1,j);
                    if c>0
                        cplot=c/maxValue;
                        lh.Color=[cplot,0,0,1];
                    else
                        cplot=c/minValue;
                        lh.Color=[0,0,cplot,1];
                    end
                    plot(x(j),y(j),'.',lh,'markersize', 15)
                    hold on
                end
            end            
        end
        function ploteyy(Mesh)
            maxValue=0;
            minValue=0;
            for i=1:size(Mesh,2)
               if max(Mesh(i).E(2,:))> maxValue
                   maxValue=max(Mesh(i).E(2,:))+eps;
               end
               if min(Mesh(i).E(2,:))< minValue
                   minValue=min(Mesh(i).E(2,:))-eps;
               end               
            end
            for i=1:size(Mesh,2)
                x=Mesh(i).x;
                y=Mesh(i).y;       
                for j=1:4
                    c=Mesh(i).E(2,j);
                    if c>0
                        cplot=c/maxValue;
                        lh.Color=[cplot,0,0,1];
                    else
                        cplot=c/minValue;
                        lh.Color=[0,0,cplot,1];
                    end
                    plot(x(j),y(j),'.',lh,'markersize', 15)
                    hold on
                end
            end            
        end  
                function plotexy(Mesh)
            maxValue=0;
            minValue=0;
            for i=1:size(Mesh,2)
               if max(Mesh(i).E(3,:))> maxValue
                   maxValue=max(Mesh(i).E(3,:))+eps;
               end
               if min(Mesh(i).E(3,:))< minValue
                   minValue=min(Mesh(i).E(3,:))-eps;
               end               
            end
            for i=1:size(Mesh,2)
                x=Mesh(i).x;
                y=Mesh(i).y;       
                for j=1:4
                    c=Mesh(i).E(3,j);
                    if c>0
                        cplot=c/maxValue;
                        lh.Color=[cplot,0,0,1];
                    else
                        cplot=c/minValue;
                        lh.Color=[0,0,cplot,1];
                    end
                    plot(x(j),y(j),'.',lh,'markersize', 30)
                    hold on
                end
            end            
        end 
        function U=getResults(X,Y,Mesh)
            flag=false;
            if isempty(Mesh(1).u)
                errordlg('Results not Parsed Yet')
            end
            for i=1:size(Mesh,2)
                x=Mesh(i).x; y=Mesh(i).y;
               if sum(X>=x)>0 && sum(X<=x)>0
                   if sum(Y>=y)>0 && sum(Y<=y)>0
                       U=Mesh(i).getU(X,Y);
                       flag=true;
                   end
               end     
             
            end
            if flag==false
                errordlg('Point Specified is outside of Mesh')  
            end
        end
        function [Z]=buildSurface(X,Y,NN,resultColumn)
            % Method Builds a surface from X and Y meshgrids a a return
            % nodal array
            Z=zeros(size(X,1),size(X,2));
            for i =1: size(X,1)
                for j=1:size(X,2)
                    for k=1:size(NN,1)
                        if (NN(k,2)==X(i,j) && NN(k,3)==Y(i,j))
                            Z(i,j)=NN(k,resultColumn);
                        end
                    end
                end
            end
        end        
        
    end
    
end

