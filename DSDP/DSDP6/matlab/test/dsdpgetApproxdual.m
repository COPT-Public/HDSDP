function [y, fub] = dsdpgetApproxdual(Amat, C)

m = length(Amat);
tol = 1e-04;

y0 = zeros(m, 1);
[emax, ~] = eigrad(Amat, C, y0);
fbest = emax;

total = 0;
baserate = 1.1;
lambda = 0.0;
maxiter = max(20, min(m / 10, 100));
baseiter = max(maxiter / 10, 20);

for i = 1:1000
    
    [y0, fub, iter] = dsdpAPLsolve(y0, Amat, C, baseiter, baserate * 2 * i / (2 * i + 1), lambda);
    total = total + iter;
    
    if (fub + lambda) > (fbest + lambda) * 2 * i / (2 * i + 1)
        baseiter = ceil(baseiter * 1.1);
        baserate = min(0.9,  baserate * 1.1);
    elseif (fub + lambda) < (fbest + lambda) * 2 * i / (2 * i + 1) * 0.8
        baserate = baserate * 0.9;
        
        if lambda > 0
            lambda = lambda * 0.9;
        else
            lambda = lambda / 0.9;
        end % End if
        
    end % End if
    
    if fub < tol
        lambda = 0.0;
    end % End if
    
    fbest = min(fub, fbest);
    
    if (fub + lambda) < tol || total >= maxiter
        break;
    end % End if
    
    % fprintf("%2d %10.3e %10.3e\n", i, fub, -lambda);
    
end % End for
y = y0;
end % End function