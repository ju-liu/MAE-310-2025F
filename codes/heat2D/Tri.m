function val = Tri(aa, xi, eta)

if aa == 1
  val = xi;
elseif aa == 2
  val = eta;
elseif aa == 3
  val = 1.0 - xi - eta;
else
  error('Error: value of a should be 1,2, or 3.');
end