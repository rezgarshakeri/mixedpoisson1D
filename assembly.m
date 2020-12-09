function [K, F] = assembly(msh, k, num_quadr_pts_in_1d, LM, g)

neq = max(max(LM));
num_elem = msh.num_elem;
% update dof per element, u+p
neldof = size(LM,1);

% the size is the total number of global dof
K = zeros(neq,neq);
Fg = zeros(neq,1);
Fb = zeros(neq,1);

for e=1:num_elem
     [ke_uu,ke_up, ke_pu, ke_pp, fe_u, fe_p] = poisson1Delem(msh, k, num_quadr_pts_in_1d, e);
     Ke = [ke_uu ke_up;...
           ke_pu ke_pp];
     Fe = [fe_u;fe_p];
    temp=LM(:,e);
    for i=1:neldof
        I=temp(i);
        if I>0
            Fb(I) = Fb(I) + Fe(i);
            for j=1:neldof
                J=temp(j);
                if J>0
                    K(I,J)=K(I,J)+ Ke(i,j);
				else
					Fg(I)=Fg(I) - (Ke(i,j)*g(j,e));
                end               
            end
        end
    end
end
F=Fg+Fb;

end

