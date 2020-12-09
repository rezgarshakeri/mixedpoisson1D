function msh = get_mesh(a, b, polydof_u, polydof_p, nelx, plot_mesh, discontinuous)

num_elem  = nelx;         % number of elements
node_dof_u = 1;         % degrees-of-freedom per node

if polydof_u == 1
% mesh specifications
el_node  = polydof_u+1;         % number of nodes per element 
lpx = nelx+1;             % number of nodes in the x direction 
% generate the IEN connectivity array
IEN = zeros(el_node, num_elem);
for i=1:num_elem
    IEN(1,i)=i;
    IEN(2,i)=i+1;
end 

elseif polydof_u == 2
% mesh specifications
el_node  = polydof_u+1;         % number of nodes per element 
lpx = 2*nelx+1;             % number of nodes in the x direction 
% generate the IEN connectivity array
IEN = zeros(el_node, num_elem);
for i=1:num_elem
    IEN(1,i)=i + (i-1);
    IEN(2,i)=i + (i-1) + 1;
    IEN(3,i)=i + (i-1) + 2;
end 

else
    error('p > 2 is not supported')
end


xu = linspace(a, b,lpx)';      % equal bisection of the x nodes
    
% get the coords of pressure nodes
if polydof_u == 1 && polydof_p == 0
   
   i=1:size(xu)-1;
   xp = (xu(i)+xu(i+1))/2;
elseif polydof_u == 1 && polydof_p == 1

   xp = xu;    

elseif polydof_u == 2 && polydof_p == 0
   
   i=2:2:size(xu)-1;
   xp = xu(i);
elseif polydof_u == 2 && polydof_p == 1
    
   i=1:2:size(xu);
   xp = xu(i);  
end

if polydof_u == 2 && discontinuous == 1
       i=1:size(xu)-1;
       xp = (xu(i)+xu(i+1))/2; 
end

nd = 0;     % number of essential boundary conditions for u
                   
ndof_u  = lpx*node_dof_u - nd;  % number of unknown 
neldof_u = el_node*node_dof_u;   % dofs per element

msh=struct();

msh.conn = IEN;
msh.coords = xu; 
msh.coords_p = xp; 
msh.num_elem = num_elem;
msh.node_dof = node_dof_u;
msh.el_node = el_node;
msh.num_nodes = lpx;
msh.ndof_u = ndof_u;
% dof of velocity per element
msh.neldof_u = neldof_u;
% dof of pressure per element
msh.neldof_p = polydof_p+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for plot

if (strcmp(plot_mesh, 'yes'))

figure(1)

    for i = 1:num_elem
        if polydof_u == 1
        XX = [xu(IEN(1,i)) xu(IEN(2,i))];
        YY = [0*i 0*i];
        plot(XX,YY, 'LineWidth',2);hold on;

        XXnode = [xu(IEN(1,i)) xu(IEN(2,i))];
        YYnode = [0*i 0*i];
        text(XXnode(1),YYnode(1),sprintf('%0.5g',IEN(1,i)),'FontSize', 16);
        text(XXnode(2),YYnode(2),sprintf('%0.5g',IEN(2,i)),'FontSize', 16);
        elseif polydof_u == 2
        XX = [xu(IEN(1,i)) xu(IEN(2,i)) xu(IEN(3,i))];
        YY = [0*i 0*i 0*i];
        plot(XX,YY, 'LineWidth',2);hold on;

        XXnode = [xu(IEN(1,i)) xu(IEN(2,i)) xu(IEN(3,i))];
        YYnode = [0*i 0*i 0*i];
        text(XXnode(1),YYnode(1),sprintf('%0.5g',IEN(1,i)),'FontSize', 16);
        text(XXnode(2),YYnode(2),sprintf('%0.5g',IEN(2,i)),'FontSize', 16);
        text(XXnode(3),YYnode(3),sprintf('%0.5g',IEN(3,i)),'FontSize', 16);

        end

    end
    title(strcat('Q', num2str(polydof_u),'P',num2str(polydof_p)))
    set(gca,'FontName','Helvetica','FontSize',16)
end

hold off


end