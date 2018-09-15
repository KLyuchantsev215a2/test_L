 function [rho] = ComputeRho(i,x,m,N,h)

rho=0;
for j = 1:N
    rho=rho+m*ComputeW(i,j,x,h); 
end

