clear; close all;

n = 1000;
tol = 1e-10;
x0 = zeros(n, 1);
ngood = 0;

rng(10);

A = sprandsym(n, 0.001);
A = A' * A + speye(n) * 1;
b = randn(n, 1) * 100;
P = inv(A);
% P = sparse(diag(diag(A)));

xtmp = (P') \ (P \ b);
tic;
[x] = dsdpConjGrad(A, b, P, x0, true);
toc;
tic;
[x] = dsdpConjGrad(A, b, P, x0, false);
toc;