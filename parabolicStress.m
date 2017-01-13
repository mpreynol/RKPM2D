function [tau] = parabolicStress(y,y0,V,h)
%Function returns the stress at point 'y' for a given beam depth and Total
%Shear load. Thickness assumed unity. Mathew Reynolds Jan 11, 2017
% y0 is the position of the neutral axis (used for coordinate trans).
tau=zeros(2,1);
if nargin==1
   y0=0.5;
   V=-3.14E1;
   h=1;
end
y=y-y0;
tau(2)=3*V/(2*h)*(1-4*y.^2/(h^2));
tau(1)=0;
% Wolfram Verification: int(3*V/(2*h)*(1-4*y^2/(h^2)),y=-h/2..h/2)
end

