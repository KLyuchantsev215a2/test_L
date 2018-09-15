function [stress2x2] = ComputeStress(F2x2,mu,k)
% Subroutine which computes Cauchy stress tensor as a function
% of the deformation gradient F
% input:  F = deformation gradient
%         mu = shear modulus
%          k = bulk modulus
% output: stress = Cauchy stresss (true stress) 
%
%   neo-Hookean material
%
F = F2x2;
F(3,3) = 1;
J = det(F);
B = F*F';  % left Cauchy-Green tensor
Biso = J^(-2/3)*B;    % isochoric part of the left Cauchy-Green tensor
devBiso = Biso - 1/3*trace(Biso)*eye(3);   % deviatoric part Biso
stress = (2*mu*devBiso + k/10*(J^5-J^(-5))*eye(3) )/J;
stress2x2(1:2,1:2) = stress(1:2,1:2);