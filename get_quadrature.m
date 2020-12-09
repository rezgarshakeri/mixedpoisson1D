
function [w,q] = get_quadrature(ngp,quadmethod)

if (strcmp(quadmethod, 'GAUSS'))
    [w, q] = GaussQuadrature(ngp);
elseif (strcmp(quadmethod, 'LGL'))
    [w, q] = LobattoQuadrature(ngp); 
end

end
  


function [W, qref1d] = GaussQuadrature(Q)
%input: Q: Number of quadrature points (Gauss)
%output:qref1d: Gauss quadrature points
%            W: Gauss weights
% Golub-Welsch algorithm: (Brute force version by Trefethen-Bau)
% to calculate Gauss points and weights using Legendre weight function 
%
    beta = 0.5./sqrt(1-(2*(1:Q-1)).^(-2));
    [V,D]=eig(diag(beta,1)+diag(beta,-1));
    [x,i]=sort(diag(D)); 
    w=2*V(1,i).^2';

    W = w';
    qref1d = x';
end

function [W, qref1d] = LobattoQuadrature(Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_Q(x). Useful for numerical
% integration and spectral methods.
%
% Reference on LGL nodes and weights:
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncation
    Q1=Q;

    % Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    x=-cos(pi*(0:Q-1)/(Q-1))';

    % The Legendre Vandermonde Matrix
    P=zeros(Q1,Q1);

    % Compute P_(Q) using the recursion relation
    % Compute its first and second derivatives and
    % update x using the Newton-Raphson method.

    xold=2;

    while max(abs(x-xold))>eps

        xold=x;

        P(:,1)=1;    P(:,2)=x;

        for k=2:Q-1
            P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
        end

        x=xold-( x.*P(:,Q1)-P(:,Q-1) )./( Q1*P(:,Q1) );

    end

    w=2./((Q-1)*Q1*P(:,Q1).^2);
    W = w';
    qref1d = x';
    
 end