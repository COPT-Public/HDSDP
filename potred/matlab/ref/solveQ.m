function [alpha] = solveQ(Q, d)

alpha = zeros(2, 1);
a = Q(1, 1); b = Q(1, 2); c = Q(2, 2);

denom = a * c - b * b;
alpha(1) = (d(1) * c - d(2) * b) / denom;
alpha(2) = (d(2) * a - d(1) * b) / denom;

end % End function