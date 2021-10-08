clear;
clc;
% load("zipcode_24_pfail_3.mat");
load("kappa_10_pfail_3.mat");
A = data.A;
b = data.b;
[m, n] = size(A);
tol = data.bestloss * 1.5;
init_x = randn(n, 1);
maxiter = 400;
batchsize = 32;
gamma = sqrt(maxiter * m / batchsize);
alpha_0 = 10;
beta = 0;
rng(123);

[sol, info] = proxpt(A, b, gamma, 0, init_x, maxiter, tol, ...
    true, 0, alpha_0, 0, true);

% [sol, info] = proxlinbatch(A, b, gamma, beta, init_x, maxiter, tol, ...
%     true, batchsize, 0, alpha_0, true);

[sol, info] = proxptbatch(A, b, gamma, init_x, maxiter, tol, ...
    true, batchsize, 0, alpha_0, true, 0);