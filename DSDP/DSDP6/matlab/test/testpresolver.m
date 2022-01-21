clear;
clc;

m = 10;
n = 20;
sps = 0.1;

rng(24);

A = cell(m, 1);

for i = 1:m
    A{i} = sprandsym(n, sps);
end % End for

b = rand(m, 1);
C = sparse(1:n, 1:n, rand(n, 1));

[pA, pb, pC, pscaler, dscaler] = dsdpPresolver(A, b, C);

cvx_begin sdp
cvx_solver mosek
variable X(n, n) symmetric
dual variable S
minimize trace(C * X)

for i = 1:m
    trace(A{i} * X) == b(i);
end % End for

X >= 0 : S;
cvx_end


cvx_begin sdp
cvx_solver mosek
variable pX(n, n) symmetric
dual variable pS
minimize trace(pC * pX)

for i = 1:m
    trace(pA{i} * pX) == pb(i);
end % End for

pX >= 0 : pS;
cvx_end

[pX, py] = dsdpPostsolver(X, y, pscaler, dscaler);