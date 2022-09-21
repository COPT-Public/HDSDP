function [x, r] = mycg(M, b, x0, tol)
% Solve M \ b via conjugate gradient

x = x0;
r = b - M * x;
d = r;
Md = M * d;
Pinvr = r;

for i = 1:100
    
    rTPinvr = r' * Pinvr;
    alpha = rTPinvr / (d' * Md);
    x = x + alpha * d;
    rnew = r - alpha * Md;
    Pinvr = rnew;
    beta = rnew' * Pinvr / rTPinvr;
    d = Pinvr + beta * d;
    Md = M * d;
    r = rnew;
    
    if norm(r) < tol
        break;
    end % End if 
    
    fprintf("%3d   %10.3e %10.3e \n", i, alpha, norm(r));
    
end % End for

fprintf("%3d %10.3e \n", i, norm(r));

end % End function