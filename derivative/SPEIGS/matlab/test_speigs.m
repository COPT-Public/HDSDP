% Test SPEIGS package
clear; clc; close all;

rng(100);
n = 2000;
opts.gthresh = 0.7;
opts.tol = 1e-10;
opts.quiet = false;

% Zero
A = sparse(n, n);
[V, e] = speigs(A, opts);
fprintf("Error %e \n", norm(V * diag(e) * V' - A, 'fro'));

% Diagonal
a = randn(n, 1);
a(a < 0) = 0;
A = sparse(diag(a));
[V, e] = speigs(A, opts);
fprintf("Error %e \n", norm(V * diag(e) * V' - A, 'fro'));

% Two-two
A = sprandsym(n, 0.00001);
[V, e] = speigs(tril(A), opts);
fprintf("Error %e \n", norm(V * diag(e) * V' - A, 'fro'));

% Rank-one
A = sparse(a * a');
[V, e] = speigs(tril(A), opts);
fprintf("Error %e \n", norm(V * diag(e) * V' - A, 'fro'));

% Submatrix
A = sprandsym(n, 0.1 / n);
Alow = tril(A);
[V, e] = speigs(Alow, opts);
fprintf("Error %e \n", norm(V * diag(e) * V' - A, 'fro'));

% General
A = sprandsym(n, 0.9);
Alow = tril(A);
[V, e] = speigs(Alow, opts);
fprintf("Error %e \n", norm(V * diag(e) * V' - A, 'fro'));

