function [penalty, infeas, ax] = dsdpgetpenaltyObj(A, X, b, bound)

m = length(b);
Ax = zeros(m, 1);

for i = 1:m
    Ax(i) = trace(A{i} * X);
end % End for

ax = Ax;
infeas = abs(Ax - b);
penalty = norm(infeas, 1) * bound;

end % End function