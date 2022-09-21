clear; clc; close all;

n = 1000;
tol = 1e-10;
x0 = zeros(n, 1);
ngood = 0;
for i = 1:1

A = sprandsym(n, 0.6) + speye(n) * 1000;
b = randn(n, 1) * 100;

[x1, e1] = drsomcg(A, b, x0, tol);
% [x1, e] = drsomcg(A, b, x0);
% [x2, e2] = drsomcg3(A, b, x0, tol);
[x3, e3] = mycg(A, b, x0, tol);

if norm(e1) <= norm(e3)
    ngood = ngood + 1;
end % End if


end % End for
