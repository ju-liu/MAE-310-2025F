clear all; clc;

kappa = 1.0;

% Exact solution
exact   = @(x,y) x*x*(1-x)*y*(1-y);
exact_x = @(x,y)  (2*x-3*x*x)*y*(1-y);
exact_y = @(x,y) x*x*(1-x)*(1-2*y);

f = @(x,y) kappa*( (6*x-2)*(y-y*y) + 2*x*x*(1-x) ); % source term

n_int = 3;
[qp, weight] = Gauss_tri( n_int );
xi = qp(:,1);
eta = qp(:,2);

mesh = read_gmsh_mesh("../../gmsh-files/square_100.m");

IEN  = mesh.tri;
x_coor = mesh.coords(:,1);
y_coor = mesh.coords(:,2);

n_el = size(IEN, 1);
n_en = size(IEN, 2);
n_np = length(x_coor);

% ID array
ID = zeros(n_np, 1);
dir_nodes = unique(mesh.edge_dirichlet(:))';
free_nodes = setdiff(1:n_np, dir_nodes);
n_eq= length(free_nodes);
ID(free_nodes) = 1 : n_eq;

LM = ID(IEN);

% sparse matrix allocation
K = spalloc(n_eq, n_eq, 9*n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );

  k_ele = zeros(n_en, n_en); % element stiffness matrix
  f_ele = zeros(n_en, 1);    % element load vector

  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Tri(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Tri(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Tri_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end

    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

    for aa = 1 : n_en
      Na = Tri(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Tri_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

      f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;

      for bb = 1 : n_en
        Nb = Tri(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Tri_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

        k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end % end of bb loop
    end % end of aa loop
  end

  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);

      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end
    end
  end
end

sol = K \ F;

disp = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    disp(ii) = sol(index);
  else
    % Dirichlet node g = 0 (nothing to do)
  end
end

% visualization
figure; hold on;

trisurf(IEN, x_coor, y_coor, disp);

shading interp;          % smooth color interpolation
colormap jet;            % color map
colorbar;                % show scale
axis equal;              % preserve geometry
view(2);                 % 2D view; comment out for 3D
title('FEM Solution');
xlabel('x'); ylabel('y');


% error calculation
n_int = 19;
[qp, weight] = Gauss_tri(n_int);
xi  = qp(:,1);
eta = qp(:,2);

L2_error  = 0.0;
H1_error  = 0.0;

for ee = 1:n_el
  % node indices of element ee
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );

  dx_dxi = 0; dx_deta = 0;
  dy_dxi = 0; dy_deta = 0;

  % gradients of shape functions (reference triangle)
  for a = 1:3
    [Na_xi, Na_eta] = Tri_grad(a, 0, 0);
    dx_dxi  = dx_dxi  + x_ele(a)*Na_xi;
    dx_deta = dx_deta + x_ele(a)*Na_eta;
    dy_dxi  = dy_dxi  + y_ele(a)*Na_xi;
    dy_deta = dy_deta + y_ele(a)*Na_eta;
  end

  J = dx_dxi*dy_deta - dx_deta*dy_dxi;

  % L2, H1 integration over triangle
  for ll = 1:n_int

    % map (xi,eta) → physical coordinates
    x_l = 0; y_l = 0;
    uh  = 0;
    uh_x = 0; uh_y = 0;

    % compute shape functions and gradients
    for a = 1:3
      Na = Tri(a, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Tri_grad(a, xi(ll), eta(ll));

      % physical coordinates
      x_l = x_l + x_ele(a)*Na;
      y_l = y_l + y_ele(a)*Na;

      % numerical solution uh
      uh = uh + disp(nodes(a))*Na;

      % physical gradients of shape function
      Na_x = ( Na_xi*dy_deta - Na_eta*dy_dxi ) / J;
      Na_y = (-Na_xi*dx_deta + Na_eta*dx_dxi ) / J;

      uh_x = uh_x + disp(nodes(a))*Na_x;
      uh_y = uh_y + disp(nodes(a))*Na_y;
    end

    % exact values
    u  = exact(x_l, y_l);
    ux = exact_x(x_l, y_l);
    uy = exact_y(x_l, y_l);

    % accumulate error
    L2_error = L2_error + weight(ll)*J * (uh - u)^2;
    H1_error = H1_error + weight(ll)*J * ((uh_x - ux)^2 + (uh_y - uy)^2);
  end
end

L2_error = sqrt(L2_error);
H1_error = sqrt(H1_error);

fprintf('L2 error  = %.6e\n', L2_error);
fprintf('H1 error  = %.6e\n', H1_error);

% EOF