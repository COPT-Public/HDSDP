clear;
% clc;
close all;

fname = fullfile('data', 'p_ADLITTLE.SIF.mps');
data = preprocess(fname);

rng(24);

A = data.A;
b = data.b;
c = data.c;

savepotdata(data);

linesearch = false;
neweigs = false;

model.A = A;
model.sense = '=';

[m, n] = size(A);

nrmb = norm(b, 1); nrmc = norm(c, 1);
b = b / (nrmb + 1); c = c / (nrmc + 1);
model.rhs = b;
model.obj = c;

grbsol = gurobi(model);

HSDAA = [sparse(m, m), A, sparse(m, n), sparse(m, 1), -b;
         -A',         sparse(n, n),  -speye(n), sparse(n, 1), c;
         b',          -c',  sparse(1, n), -1, 0];
     
% HSDAA = HSDAA' * inv(HSDAA * HSDAA') * HSDAA; 
[D, E, HSDA] = ruizscale(HSDAA, 30);
lpsol = potreduceLp(HSDA, m, 100, false, linesearch, neweigs, 1);
sol = lpsol .* E;

kappa = sol(end - 1);
tau = sol(end);
y = sol(1:m);
y = y / tau;
s = sol(m + n + 1 : m + 2 * n) / tau;
x = sol(m + 1: m + n) / tau;

fprintf("%20s %10.3e, %10.3e  %10.3e  %10.3e \n", fname, c' * x, ...
    norm(A * x - b) / (1 + norm(b, 1)), norm(A' * y + s - c) / (norm(c, 1) + 1), ...
    (c' * x - b' * y) / (abs(c' * x) + abs(b' * y) + 1));
