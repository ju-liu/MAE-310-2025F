clear all; clc;

kappa = 1.0;

% Exact solution
exact   = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y)  (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

num_int = 3;
[qp, qw] = Gauss_tri(num_int);

mesh = read_gmsh_mesh("../../../square.m");

IEN  = mesh.tri;
coor_x = mesh.coords(:,1);
coor_y = mesh.coords(:,2);

n_el = size(IEN, 1);
n_en = size(IEN, 2);
n_np = size(coor_x, 1);

% ID array
ID = zeros(n_np, 1);
dir_nodes = unique(mesh.edge_dirichlet(:))';
free_nodes = setdiff(1:n_np, dir_nodes);
ID(free_nodes) = 1 : length(free_nodes);

% EOF