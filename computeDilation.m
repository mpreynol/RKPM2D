function [aArray] = computeDilation(NN)
    % Matrix with each element a data structure with distance vector between points i and j (sym).
    numberOfNodes=size(NN,1);
    A=zeros(numberOfNodes,numberOfNodes,3);
    aArray=ones(obj.numberOfNodes,1)*inf;
    for i=1:numberOfNodes
        for j=1:numberOfNodes
            if i~=j
                A(i,j,1)=NN(i,2)-NN(j,2);
                A(i,j,2)=NN(i,3)-NN(j,3);
                A(i,j,3)=sqrt(A(i,j,1)^2+A(i,j,2)^2);
            end
        end
    end
    for i=1:numberOfNodes
        for j=1:numberOfNodes
            if i~=j
                if A(i,j,1)~=0 && A(i,j,2)~=0 % Point is not CoPlanar
                    if A(i,j,3)<aArray(i)
                        aArray(i)=A(i,j,3);
                    end
                end
            end
        end
    end
end


