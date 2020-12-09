function [ID, LM, g] = get_ID_LM(msh, discontinuous)


num_nodes = msh.num_nodes;
el_node = msh.el_node;
num_elem = msh.num_elem;
% we update dof per element in assembly.m
neldof_u = msh.neldof_u;
neldof_p = msh.neldof_p;
ndof_u = msh.ndof_u;

IEN = msh.conn;

ID = 1:num_nodes;

xu = msh.coords;
xp = msh.coords_p;

% u(0) = 0, u(1) = sin(1)
e_bc=zeros(size(ID));
e_bc(1) = 0;
e_bc(end) = 0;

% before adding pressure dof
LM1 = zeros(el_node,num_elem);
for i = 1:num_elem
    for j = 1:el_node
        idd = ID(IEN(j,i));
        LM1(j,i) = idd;
    end
end

g1 = zeros(el_node,num_elem);
for i = 1:num_elem
    for j = 1:el_node
        idd1 = e_bc(IEN(j,i));
        g1(j,i) = idd1;
    end
end

% this is for neldof_p = 2, continuous case to update LM1
IENp = zeros(2, num_elem);
for i=1:num_elem
    IENp(1,i)=i;
    IENp(2,i)=i+1;
end 

% updating LM after adding pressure nodes (nelpdof)
ndof = ndof_u;
if neldof_p == 1 
    LM = LM1;
    g = g1;
    for i = 1 :num_elem
        for j = 1:neldof_p
    LM(neldof_u+j,i) = ndof + j;
    g(neldof_u+j,i) = 0;
        end
        ndof = ndof+j;
    end

  elseif  neldof_p == 2 && discontinuous==1
    LM = LM1;
    g = g1;
    for i = 1 :num_elem
        for j = 1:neldof_p
    LM(neldof_u+j,i) = ndof + j;
    g(neldof_u+j,i) = 0;
        end
        ndof = ndof+j;
    end
    
 elseif  neldof_p == 2 && discontinuous==0
    LM = LM1;
    g = g1;
    for i = 1 :num_elem
        for j = 1:neldof_p
    LM(neldof_u+j,i) = ndof_u + IENp(j,i);
    g(neldof_u+j,i) = 0;
        end
    end
    
end

end