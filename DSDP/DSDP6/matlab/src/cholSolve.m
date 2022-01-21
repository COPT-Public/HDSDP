function [X] = cholSolve(L, p, pinv, B)
% Solve linear system A * X = B, A = L * L'

B = B(p, :);
X1 = L \ B;
X = (L') \ X1;
X = X(pinv, :);

end % End function