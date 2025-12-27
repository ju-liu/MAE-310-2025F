clear all; clc;

% define material parameters
% kappa_ij = kappa delta_ij
kappa = 1.0;

% exact solution
exact   = @(x,y) sin(x)*sin(y);
exact_x = @(x,y) cos(x)*sin(y);
exact_y = @(x,y) sin(x)*cos(y);
g_data  = @(x,y) exact(x,y);
h_data  = @(x,y) -kappa * cos(x)*sin(y);

f = @(x,y) kappa * ( 2.0*sin(x)*sin(y) );

% quadrature rule
n_int = 3;
[xi, eta, weight] = Gauss_tri(n_int);

% load mesh
mesh = read_gmsh_mesh('../../gmsh-files/square_50.m');

IEN    = mesh.tri;
x_coor = mesh.coords(:,1);
y_coor = mesh.coords(:,2);

n_el = size(IEN, 1);
n_en = size(IEN, 2);
n_np = length(x_coor);

dir_nodes  = unique(mesh.edge_dirichlet(:))';
free_nodes = setdiff(1:n_np, dir_nodes);

n_eq = length(free_nodes);

% ID, and LM
ID = zeros(n_np, 1);
ID(free_nodes) = 1 : n_eq;
LM = ID(IEN);

% sparse matrix allocation
K = spalloc(n_eq, n_eq, 9*n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
  x_ele = x_coor(IEN(ee,:));
  y_ele = y_coor(IEN(ee,:));

  k_ele = zeros(n_en, n_en);
  f_ele = zeros(n_en, 1);

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

    if detJ <= 0
      error('Non-positive Jacobian at element %d, qp %d', ee, ll);
    end

    for aa = 1 : n_en
      Na = Tri(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Tri_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      f_ele(aa) = f_ele(aa) + weight(ll) * detJ * Na * f(x_l, y_l);
      for bb = 1 : n_en
        [Nb_xi, Nb_eta] = Tri_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

        k_ele(aa,bb) = k_ele(aa,bb) + weight(ll)* detJ * ...
          kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end
    end
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
          nodeQ = IEN(ee, bb);
          F(PP) = F(PP) - k_ele(aa,bb) * g_data(x_coor(nodeQ), y_coor(nodeQ));
        end
      end
    end
  end
end

% Neumann boundary condition
IEN_h   = mesh.edge_neumann;
n_el_h  = size(IEN_h, 1);
n_en_h  = 2;
n_int_h = 3;
[ss,ww] = Gauss_1D(n_int_h, -1, 1);

for ee = 1 : n_el_h
  x_ele = x_coor(IEN_h(ee,:));
  y_ele = y_coor(IEN_h(ee,:));

  for ll = 1 : n_int_h
    x_l = 0.0; y_l = 0.0;
    dx_ds = 0.0; dy_ds = 0.0;
    for aa = 1 : n_en_h
      x_l = x_l + PolyShape(1, aa, ss(ll), 0) * x_ele(aa);
      y_l = y_l + PolyShape(1, aa, ss(ll), 0) * y_ele(aa);

      dx_ds = dx_ds + PolyShape(1, aa, ss(ll), 1) * x_ele(aa);
      dy_ds = dy_ds + PolyShape(1, aa, ss(ll), 1) * y_ele(aa);
    end

    detJ = sqrt(dx_ds^2 + dy_ds^2);

    hh = h_data(x_l, y_l);

    for aa = 1 : n_en_h
      PP = ID(IEN_h(ee,aa));
      if PP > 0
        F(PP) = F(PP) + ww(ll) * detJ * PolyShape(1, aa, ss(ll), 0) * hh;
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
    disp(ii) = g_data(x_coor(ii), y_coor(ii));
  end
end

% visualize
trisurf(IEN, x_coor, y_coor, disp);
shading interp;
colormap jet;
axis equal;
view(2);

% error calculator
L2_error = 0.0;
H1_error = 0.0;

n_int = 19;
[xi, eta, weight] = Gauss_tri(n_int);

for ee = 1 : n_el
  x_ele = x_coor(IEN(ee,:));
  y_ele = y_coor(IEN(ee,:));
  u_ele = disp(IEN(ee,:));

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

    uh = 0.0; uh_x = 0.0; uh_y = 0.0;
    for aa = 1 : n_en
      Na = Tri(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Tri_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

      uh   = uh   + u_ele(aa) * Na;
      uh_x = uh_x + u_ele(aa) * Na_x;
      uh_y = uh_y + u_ele(aa) * Na_y;
    end

    u  = exact(x_l, y_l);
    ux = exact_x(x_l, y_l);
    uy = exact_y(x_l, y_l);

    L2_error = L2_error + weight(ll) * detJ * (uh - u)^2;
    H1_error = H1_error + weight(ll) * detJ *...
      ( (uh_x - ux)^2 + (uh_y - uy)^2 );
  end
end

L2_error = sqrt(L2_error);
H1_error = sqrt(H1_error);

fprintf('L2 error = %.6e\n', L2_error);
fprintf('H1 error = %.6e\n', H1_error);

% EOF