function [emin, grad] = eigrad(A, C, y)
% Evaluate gradient of lambda_min for sub-gradient method
m = length(y);

S = C - dsdpgetATy(A, y);
[v, emin] = eigs(S, 1, 'smallestreal');

grad = zeros(m, 1);
if emin > 0
    return;
else
    for i = 1:m
        grad(i) = - v' * A{i} * v;
    end % End for
end % End if 

end % End function