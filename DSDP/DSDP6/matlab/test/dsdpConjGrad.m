function [x] = dsdpConjGrad(M, b, P, x0)
% Solve M \ b via conjugate gradient

x = x0;
r = b - M * x;
d = P \ r;
Md = M * d;
Pinvr = P \ r;

for i = 1:100
    rTPinvr = r' * Pinvr;
    alpha = rTPinvr / (d' * Md);
    x = x + alpha * d;
    rnew = r - alpha * Md;
    Pinvr = P \ rnew;
    beta = rnew' * Pinvr / rTPinvr;
    d = Pinvr + beta * d;
    Md = M * d;
    r = rnew;
    
    if (norm(r) < 1e-06)
        fprintf("%d %10.3e \n", i, norm(r));
        break;
    end % End if 
    
    fprintf("%d %10.6e %20.10e \n", i, alpha, norm(r));
    
end % End for


% End function