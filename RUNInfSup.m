
clear all;
close all; 

a = 0;
b = 1;
k = 1;

% dof of velocity shape function, can be 1 or 2
polydof_u = 1;
% dof of pressure shape function, can be 0 (discontinuous) or 1 (dis/continuous pressure)
polydof_p = 0;
% to make discontinuous this must be ON otherwise must be OFF
discontinuous = 'ON';
% quadrature points can be GAUSS or LGL
quadmethod = 'GAUSS';
% num of quadrature points
num_quadr_pts = 3;

plot_mesh = 'no'; 

% loop over elements to plot convergence rate
for i=2:9
% number of element    
nelx = i;
msh = get_mesh(a, b, polydof_u, polydof_p, nelx, plot_mesh, discontinuous);
xu = msh.coords;
xp = msh.coords_p;

[ID, LM, g] = get_ID_LM(msh, discontinuous);

[K_div, K, F] = assembly(msh, k, num_quadr_pts, LM, g, discontinuous, quadmethod);
% 
% pressure dof start counting after ndof_u
ndof_u = msh.ndof_u;

% add flux
F(1) = F(1) + pressure(xu(1));
F(ndof_u) = F(ndof_u) - pressure(xu(end));


S = K(1:ndof_u,1:ndof_u);           % K_uu, a(u,v)=int(Kinv*u,v), Kinv = 1, here
AT = K(ndof_u+1:end,1:ndof_u);      % K_pu, b(u,q)
A = K(1:ndof_u,ndof_u+1:end);       % K_up, b(v,p)
C = K(ndof_u+1:end,ndof_u+1:end);   % K_pp, <p,q>_L2 which is zero in general, but I need it for inf-sup test
div_uv = K_div(1:ndof_u,1:ndof_u);  % this is (div(v),div(u))
% Inf-Sup constant from eq. (3.2) of Arnold 2009
% <u,v> + b(v,p) + b(u,q) = -lambda1 <p,q>
% I have created the equivalent of 3.2 as
H = S + div_uv; % since Kinv = 1, <u,v>_L2 = S=a(u,v)
LHS1 = [H A;
        AT 0*C];
RHS1 = -[0*H 0*A;
         0*AT C];
     
lambda1 = eig(LHS1,RHS1);
beta = sqrt(min(lambda1));

% Coercivity constant from eq. (3.4) of Arnold 2009
% a(u,v) + b(v,p) + b(u,q) = lambda2 <u,v>
% I have created the equivalent of 3.4 as
LHS2 = [S A;
        AT 0*C];
RHS2 = [H 0*A;
        0*AT 0*C];
     
lambda2 = eig(LHS2,RHS2);
alpha = abs(min(lambda2));

% mesh size, i start from 2
h = 1/i;
xx(i-1) = h;
% inf-sup constant
yy1(i-1) = beta;
% coercivity constant
yy2(i-1) = alpha;


end

figure(2)
loglog(xx,yy1,'LineWidth',3)
xlabel('log(h)','fontsize',18)
ylabel('log(\beta_h)','fontsize',18)
title(strcat('Inf-Sup Constant for mixedpoisson-','Q', num2str(polydof_u),'P',num2str(polydof_p)))
grid on
set(gca,'FontName','Helvetica','FontSize',18)

figure(3)
loglog(xx,yy2,'LineWidth',3)
title(strcat('Coercivity Constant for mixedpoisson-','Q', num2str(polydof_u),'P',num2str(polydof_p)))
xlabel('log(h)','fontsize',18)
ylabel('log(\alpha_h)','fontsize',18)
grid on
set(gca,'FontName','Helvetica','FontSize',18)

