function [ w ] = myfunWeight( x,xI,a )
z=abs(x-xI)/a;
if z<=1/2
    w=2/3-4*z^2+4*z^3;
elseif z<=1
    w=4/3-4*z+4*z^2-4/3*z^3;
else 
    w=0;
end
end

