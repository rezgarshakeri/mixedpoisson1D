 function [ke_div_uv,ke_uu,ke_up, ke_pu, ke_pp, fe_u, fe_p] = poisson1Delem(msh, k, num_quadr_pts, discontinuous,quadmethod, e)

IEN = msh.conn;
vtx_coords = msh.coords; 
neldof_u = msh.neldof_u;
neldof_p = msh.neldof_p;

kinv = 1/k;
ke_div_uv  = zeros(neldof_u,neldof_u); 
ke_uu  = zeros(neldof_u,neldof_u); 
ke_up  = zeros(neldof_u,neldof_p);    
ke_pu  = zeros(neldof_p, neldof_u);    
ke_pp = zeros(neldof_p, neldof_p);

fe_u = zeros(neldof_u,1);
fe_p = zeros(neldof_p,1);

% get coordinates of element nodes 
je = IEN(:,e);  
coords  = vtx_coords(je); 
[w , gp]  = get_quadrature(num_quadr_pts,quadmethod);           % quad points and weights 

for i = 1:num_quadr_pts  
    
  xi = gp(i);
  [Nu, Np  ]     = get_shape(xi, neldof_u, neldof_p, discontinuous);        % shape functions matrix
  [B, div, detJ]     = get_shape_dirv(xi, neldof_u, coords);     % derivative of the shape functions
      
  ke_uu = ke_uu + w(i)*Nu'*kinv*Nu*detJ;   
  ke_up = ke_up + w(i)*div'*(-1)*Np*detJ;
  ke_pu = ke_pu + w(i)*Np'*(-1)*div*detJ;
  ke_pp = ke_pp + w(i)*Np'*1*Np*detJ;  % I need this for inf-sup test <p,q>_L2
  ke_div_uv = ke_div_uv + w(i)*div'*div*detJ;   % I need this for inf-sup test <div_u,div_v>_L2

  xe = Nu*coords;  % for later if force is function of x
% 
  fp = forcing(xe,k);
  fu = 0;
  fe_u = fe_u + w(i)*Nu'*fu*detJ;    
  fe_p = fe_p + w(i)*Np'*(-fp)*detJ;
end

 end
 
