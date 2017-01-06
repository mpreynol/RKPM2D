xI=2;
x=[rand()*2;rand()*2];
Value=0;
for i=1:PointCloud.numberOfNodes
   Value=Value+PointCloud.Nodes(i).sF.getValue(x)*xI;
end
Value-xI;