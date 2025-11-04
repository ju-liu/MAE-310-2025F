function plotdisp(IEN, x_coor, disp, exact_handle)
hasExact = (nargin >= 4) && isa(exact_handle, 'function_handle');
[n_en, n_el] = size(IEN);
pp = n_en - 1;

n_sam = 20;
xi_sam = -1 : (2/n_sam) : 1;

x_sam = zeros(n_el * n_sam + 1, 1);
y_sam = x_sam; % store the exact solution value at sampling points
u_sam = x_sam; % store the numerical solution value at sampling pts

for ee = 1 : n_el
  x_ele = x_coor( IEN(:, ee) );
  u_ele = disp( IEN(:, ee) );

  if ee == n_el
    n_sam_end = n_sam+1;
  else
    n_sam_end = n_sam;
  end

  for ll = 1 : n_sam_end
    x_l = 0.0;
    u_l = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
      u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
    end

    x_sam( (ee-1)*n_sam + ll ) = x_l;
    u_sam( (ee-1)*n_sam + ll ) = u_l;
    if hasExact
      y_sam( (ee-1)*n_sam + ll ) = exact_handle(x_l);
    end
  end
end

clf;
p1 = plot(x_sam, u_sam, '--r', 'LineWidth', 2);
hold on;

if hasExact
  p2 = plot(x_sam, y_sam, '-k', 'LineWidth', 2);
  legend([p1 p2], {'FEM', 'Exact'}, 'Location', 'best');
else
  legend(p1, 'FEM', 'Location', 'best');
end

xlabel('x'); ylabel('u(x)'); grid on;
title(sprintf('FEM solution (p=%d, elements=%d)', pp, n_el));

end