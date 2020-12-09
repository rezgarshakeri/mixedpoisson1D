
function [w,gp] = gauss(ngp)
%input: ngp: Number of quadrature points (Gauss)
%output:qref1d: Gauss quadrature points
%       W: Gauss weights
% Golub-Welsch algorithm: (Brute force version by Trefethen-Bau)
% to calculate Gauss points and weights using Legendre weight function 
%
    beta = 0.5./sqrt(1-(2*(1:ngp-1)).^(-2));
    [V,D]=eig(diag(beta,1)+diag(beta,-1));
    [x,i]=sort(diag(D)); 
    w=2*V(1,i).^2;
    gp = x';
end
  
