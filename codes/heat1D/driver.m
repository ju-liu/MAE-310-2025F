clear all; clc;

% -----------------------------------------------
% material parameters and input da
% define the manufactured solution
exact   = @(x) sin(x);
exact_x = @(x) cos(x);

% define the input for the differential equation
f = @(x) sin(x);  % u,xx + f = 0
g = sin(1);       %     u(1) = g
h = -1;           %  -u,x(0) = h

% define the discretization parameters
n_el = 5;         % number of elements
n_en = 2;         % number of element nodes
n_np = n_el + 1;  % number of nodal points
hh   = 1 / n_el;  % uniform element length

x_coor = 0 : hh : 1; % nodal coordinates

% Data structures
IEN = zeros(n_en, n_el);

for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(aa, ee) = ee - 1 + aa;
  end
end

ID = 1 : n_np;
ID(end) = 0;


















% EOF