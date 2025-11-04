clear all; clc; % clean the memory and the screen

% -----------------------------------------------
% material parameters and input data
% define the manufactured solution
exact   = @(x) sin(x);
exact_x = @(x) cos(x);

% define the input for the differential equation
f = @(x) sin(x);        % u,xx + f = 0
g = exact(1);           %     u(1) = g
h = -exact_x(0);        %  -u,x(0) = h

% define the discretization parameters
pp    = 1;              % degree of shape function
n_el  = 5;              % number of elements
n_en  = pp+1;           % number of element nodes
n_np  = n_el * pp + 1;  % number of nodal points
n_eq  = n_np - 1;       % number of equations
n_int = 3;              % number of quadrature points
hh    = 1 / (n_np -1);  % length between two adjacent nodes
x_coor = 0 : hh : 1;    % nodal coordinates

% Quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

% Data structures
ID = 1 : n_np;
ID(end) = 0;

IEN = zeros(n_en, n_el); LM = IEN;

for ee = 1 : n_el
  IEN(:, ee) = (ee - 1)*pp + (1:n_en);
  LM(:, ee) = ID(IEN(:,ee));
end

% Matrix and vector allocation
K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
F = zeros(n_eq, 1);
for ee = 1 : n_el
  k_ele    = zeros(n_en, n_en);
  f_ele    = zeros(n_en, 1);
  x_ele(:) = x_coor( IEN(:, ee) );

  for ll = 1 : n_int
    x_l    = 0.0; dx_dxi = 0.0;
    for aa = 1 : n_en
      x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(ll), 0);
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(ll), 1);
    end
    dxi_dx = 1.0 / dx_dxi;

    for aa = 1 : n_en
      f_ele(aa) = f_ele(aa) + weight(ll) * PolyShape(pp, aa, xi(ll), 0) * f(x_l) * dx_dxi;
      for bb = 1 : n_en
        k_ele(aa, bb) = k_ele(aa, bb) + weight(ll) * PolyShape(pp, aa, xi(ll), 1) ...
          * PolyShape(pp, bb, xi(ll), 1) * dxi_dx;
      end
    end
  end % end of quadrature loop

  for aa = 1 : n_en
    P = LM(aa, ee);
    if P > 0
      F(P) = F(P) + f_ele(aa);
      for bb = 1 : n_en
        Q = LM(bb, ee);
        if Q > 0
          K(P,Q) = K(P,Q) + k_ele(aa, bb);
        else
          F(P) = F(P) - k_ele(aa, bb) * g;
        end
      end
    end
  end
end % end of element loop

F(1) = F(1) + h;

% Solve the matrix problem K d_temp = F
d_temp = K \ F;

% put Dirichlet value into the disp vector
disp = [d_temp; g];

% visulization
plotdisp(IEN, x_coor, disp, exact);


% EOF