clear; clc; close all;


m = 50;
n = 5000;

A1 = randn(m, n);
A = A1' * inv(A1 * A1') * A1;
gamma = norm(A1, 'fro');

% cvx_begin
% cvx_solver sedumi
% variable x(n, 1) nonnegative;
% minimize 0
% subject to 
%         A1 * x == zeros(m, 1);
%         sum(x) == 1;
% cvx_end

% [x] = firstordPot(A, 1);
% [xpot] = potreduce(A1, gamma, true);
[xpot] = potreduce(A1, gamma, false, true);
[xpot] = potreduce(A1, gamma, false, false);