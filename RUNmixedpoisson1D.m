
clear all;
close all; 

a = 0;
b = 1;
k = 1;
tol = 1e-13;
maxit = 1000;

% dof of velocity shape function, can be 1 or 2
polydof_u = 2;
% dof of pressure shape function, can be 0 (discontinuous) or 1 (dis/continuous pressure)
polydof_p = 1;
% to make discontinuous this must be 1 otherwise 0
discontinuous = 1;

plot_mesh = 'yes'; 

for i=2:8
nelx = i;
num_quadr_pts_in_1d = 3;
msh = get_mesh(a, b, polydof_u, polydof_p, nelx, plot_mesh, discontinuous);
xu = msh.coords;
xp = msh.coords_p;

[ID, LM, g] = get_ID_LM(msh, discontinuous);

[K, F] = assembly(msh, k, num_quadr_pts_in_1d, LM, g);
% 
% pressure dof start counting after ndof_u
ndof_u = msh.ndof_u;

% add flux
F(1) = F(1) + pressure(xu(1));
F(ndof_u) = F(ndof_u) - pressure(xu(end));


H = K(1:ndof_u,1:ndof_u);           % K_uu
AT = K(ndof_u+1:end,1:ndof_u);      % K_pu
A = K(1:ndof_u,ndof_u+1:end);       % K_up

S = AT*inv(H)*A;
lambda = eig(S);
U0 = zeros(size(F));
omega = 1/max(lambda);

[U, norm_res,norm_sol] = Uzawa(K,F, ndof_u, U0, omega, tol);

% % Uzawa's method
Uu = U(1:ndof_u,1);
Up = U(ndof_u+1:end,1);

Uh = Uu;
Ph = Up;
  
ue = velocity(xu,k);
pe = pressure(xp);

er_U = abs(ue - Uh);
er_P = abs(pe - Ph);

error_U=norm(ue - Uh);
error_P=norm(pe - Ph);

EP(i-1) = error_P;
EU(i-1) = error_U;



P = [H 0*A;0*AT S];
d = minres(K,F,tol,maxit,P);

% minres solution
du = d(1:ndof_u,1);
dp = d(ndof_u+1:end,1);
    
uh = du;
ph = dp;
  
ue = velocity(xu,k);
pe = pressure(xp);

er_u = abs(ue - uh);
er_p = abs(pe - ph);

error_u=norm(ue - uh);
error_p=norm(pe - ph);

Ep(i-1) = error_p;
Eu(i-1) = error_u;
h(i-1) = 1/i;

end

lglg_factor_1 = 0.05;
lglg_pwr_1 = h.^1;
lglg_factor_2 = 0.010;
lglg_pwr_2 = h.^2;
leg_enry_2 = 'O(h^1)';
leg_enry_3 = 'O(h^2)';

figure(2)
loglog(h,Ep,h,Eu,'LineWidth',2)
title(strcat('Q', num2str(polydof_u),'P',num2str(polydof_p)))
xlabel('h','fontsize',16)
ylabel('absolute error','fontsize',16)
grid on
%legend('Matrix Free FEM')
set(gca,'FontName','Helvetica','FontSize',14)
hold on
loglog(h,lglg_factor_1*(lglg_pwr_1),'r:',h, lglg_factor_2*(lglg_pwr_2),'b--', 'LineWidth',2);
legend('pressure','velocity',leg_enry_2,leg_enry_3,'Location','northwest')

figure(3)
plot(er_u)
title(strcat('Q', num2str(polydof_u),'P',num2str(polydof_p), '-This is the nodal value error for velocity: (ue - uh)'))
set(gca,'FontName','Helvetica','FontSize',16)

figure(4)
plot(er_p)
title(strcat('Q', num2str(polydof_u),'P',num2str(polydof_p), '-This is the nodal value error for pressure: (pe - ph)'))
    set(gca,'FontName','Helvetica','FontSize',16)

figure(5)
plot(er_U)
title(strcat('Q', num2str(polydof_u),'P',num2str(polydof_p), '-(ue - uh) with Uzawa Method'))
set(gca,'FontName','Helvetica','FontSize',16)

figure(6)
plot(er_P)
title(strcat('Q', num2str(polydof_u),'P',num2str(polydof_p), '-(pe - ph) with Uzawa Method'))
    set(gca,'FontName','Helvetica','FontSize',16)
    
    
    