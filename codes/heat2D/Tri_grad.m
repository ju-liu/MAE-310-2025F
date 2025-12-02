function [val_xi, val_eta] = Tri_grad(aa, xi, eta)

if aa == 1
  val_xi  = 1.0;
  val_eta = 0.0;
elseif aa == 2
  val_xi  = 0.0;
  val_eta = 1.0;
elseif aa == 3
  val_xi  = -1.0;
  val_eta = -1.0;
else
  error('Error: value of a should be 1,2, or 3.');
end