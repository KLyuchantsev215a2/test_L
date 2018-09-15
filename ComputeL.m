function [Ln1]=ComputeL(L,i,j,beta,m,v,rho,nabla_W)

L(1,beta,i)=L(1,beta,i)+(m/rho(1,j)*(v(1,1,j)-v(1,1,i))*nabla_W(beta));  
L(2,beta,i)=L(2,beta,i)+(m/rho(1,j)*(v(1,2,j)-v(1,2,i))*nabla_W(beta)); 

Ln1=L;