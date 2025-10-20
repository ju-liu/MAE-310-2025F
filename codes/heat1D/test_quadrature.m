clear all;
clc;
n_int = 2;
[xi, wq] = Gauss(n_int, -1, 1);

sum = 0.0;
for qua = 1 : n_int
    sum = sum + wq(qua) * xi(qua)^4;
end

syms x;
exact = int(x^4, x, -1, 1);