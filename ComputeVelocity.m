function [velocity]=ComputeVelocity(i,j,beta,dt,m,v,rho,SIG,PI,nabla_W)
   
v(1,1,i)=v(1,1,i)-dt*(m*(SIG(1,beta,i)/rho(1,i)^2+SIG(1,beta,j)/rho(1,j)^2)+PI*(1==beta))*nabla_W(beta); 
v(1,2,i)=v(1,2,i)-dt*(m*(SIG(2,beta,i)/rho(1,i)^2+SIG(2,beta,j)/rho(1,j)^2)+PI*(2==beta))*nabla_W(beta);
    
velocity=v;