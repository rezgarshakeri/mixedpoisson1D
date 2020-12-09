function [Nu, Np] = get_shape(xi, neldof_u, neldof_p)

% 1---2---3

if neldof_u == 2
    N1 = (1-xi)/2;
    N2 = (1+xi)/2;
    Nu = [N1 N2];
elseif neldof_u == 3
    N1 = xi*(xi-1)/2;
    N2 = 1-xi^2;
    N3 = xi*(xi+1)/2;
    Nu = [N1 N2 N3];
else
    error('does not support p>2');
end

% set pressure interpolation
    if (neldof_p==1)
        Np = 1;
    elseif (neldof_p==2)
        N1 = (1-xi)/2;
        N2 = (1+xi)/2;
        Np = [N1 N2];
    end

    
end
