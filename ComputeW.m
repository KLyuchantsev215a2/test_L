function [W]=ComputeW(i,j,x,h)
%kernel of SPH as a function
% of the coordinate particle i and j

% input:  %a number initial particle
          %b namber calculated particle 
          %h  blurring radius
          %x coordinate all particle 
% output: W = the force of the particle j effect on the initial i

W=0;

r=zeros(1, 2);

r(1,1)=x(1,1,i)-x(1,1,j);
r(1,2)=x(1,2,i)-x(1,2,j);

q=norm(r,2)/h;
C=1/(pi*h*h);

if (q>=0) && (q<=1)
    W=C*(15 / 7)*(2 / 3 - q*q + 1/2*q*q*q);
elseif (q >= 1) && (q <= 2)
    W = C*(5 / 14)*(2 - q)^3;
end
