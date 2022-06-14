clear;
clc;

A = sparse(10, 10);
A(1, 5) = 0.5;
A(2, 2) = 1;
A = A + A';

t = sqrt(2);
[V, E] = eigs(A);

v1 = zeros(10, 1);
v2 = zeros(10, 1);
v1(1) = t / 2;
v1(5) = t / 2;
e1 = 0.5;
v2(1) = -t / 2;
v2(5) = t / 2;
e2 = -0.5;