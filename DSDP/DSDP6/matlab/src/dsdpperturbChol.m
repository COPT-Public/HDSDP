function [R] = dsdpperturbChol(X, eps)
% Perturb the cholesky decompostion of S when necessary

if nargin < 2
    eps = 1e-06;
end % End if

[n, ~] = size(X);

if eps > 0
    
    for i = 1:10
        try
            R = lchol(X)';
            break;
        catch
            X = X + speye(n) * eps;
        end % End try
    end % End for
    
else
    R = lchol(X)';
end % End if


end % End function