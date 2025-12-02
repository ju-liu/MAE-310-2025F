clear all; clc;

% define material parameters
% kappa_ij = kappa delta_ij
kappa = 1.0;

% exact solution
exact = @(x,y) x*x*(1-x)*y*(1-y);
exact_x = @(x,y) (2*x-3*x*x)*y*(1-y);
exact_y = @(x,y) x*x*(1-x)*(1-2*y);

f = @(x,y) kappa * ( (6*x-2)*(y-y*y) + 2*x*x*(1-x) );

% quadrature rule
n_int = 3;
[xi, eta, weight] = Gauss_tri(n_int);

% load mesh
mesh = read_gmsh_mesh('../../gmsh-files/square_2.m');

IEN = mesh.tri;
x_coor = mesh.coords(:,1);
y_coor = mesh.coords(:,2);

dir_nodes = unique(mesh.edge_dirichlet(:))';
free_nodes = setdiff(1:n_np, dir_nodes);

n_el = size(IEN, 1);
n_en = size(IEN, 2);
n_np = length(x_coor);
n_eq = length(free_nodes);

% ID, and LM
ID = zeros(n_np, 1);
ID(free_nodes) = 1 : n_eq;

LM = ID(IEN);

% sparse matrix allocation
K = spalloc(n_eq, n_eq, 9);
F = zeros(n_eq, 1);

for ee = 1 : n_el
  x_ele = x_coor(IEN(ee,:));
  y_ele = y_coor(IEN(ee,:));

  k_ele = zeros(n_en, n_en);
  f_ele = zeros(n_en, 1);

  for ll = 1 : n_int
    x_l = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Tri(aa, xi(ll), eta(ll));
    end


  end


end










% EOF