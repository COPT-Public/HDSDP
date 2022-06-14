function [emax, grad] = eigrad(A, C, y)
% Evaluate gradient of lambda_min for sub-gradient method
m = length(y);
[n, ~] = size(C);

S = dsdpgetATy(A, y) - C;
[v, emax] = eigs(S, 1, 'largestreal', 'SubspaceDimension', min(40, n), 'FailureTreatment', 'keep');

grad = zeros(m, 1);

if isnan(emax)
    error("Eigenvalue computation failed");
end % End if

if emax <= -100
    return;
else
    for i = 1:m
        grad(i) = v' * A{i} * v;
    end % End for
end % End if 

end % End function