function [nabla_W]=Compute_nabla_W(i,j,x,h,beta)
%kernel of SPH as a function
% of the coordinate particle i and j

% input:  %a number initial particle
          %b namber calculated particle 
          %h  blurring radius
          %x coordinate all particle 
          %betta normal for the derivative 
% output: W = the force of the particle j effect on the initial i

nabla_W=zeros(1, 2);
    
r=zeros(1, 2);

r(1,1)=x(1,1,i)-x(1,1,j);
r(1,2)=x(1,2,i)-x(1,2,j);

q=norm(r,2)/h;
C=1/(pi*h*h);

if (q>0) && (q<1)
  nabla_W(1,beta)=C*(15 / 7)*(-2*q+3/2*q*q)*(r(1,beta)/(h*norm(r)));
elseif (q >= 1) && (q <= 2)
   nabla_W(1,beta) = C*(15 / 7)*(-1/2)*(2 - q)^2*(r(1,beta)/(h*norm(r)));
end
