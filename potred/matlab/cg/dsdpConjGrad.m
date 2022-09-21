function [x] = dsdpConjGrad(M, b, P, x0, restart)
% Solve M \ b via conjugate gradient

n = size(M, 1);
x = x0;
r = b - M * x;
d = P \ r;
Md = M * d;
Pinvr = P \ r;

for i = 1:100
    
    rTPinvr = r' * Pinvr;
    alpha = rTPinvr / (d' * Md);
    x = x + alpha * d;
    
    if mod(i, 20) == 1 && restart
        r = b - M * x;
        d = P \ r;
        Md = M * d;
        Pinvr = P \ r;
        continue;
    end % End if
    
    rnew = r - alpha * Md;
    Pinvr = P \ rnew;
    beta = rnew' * Pinvr / rTPinvr;
    d = Pinvr + beta * d;
    Md = M * d;
    
    r = rnew;
    
    if (norm(r) < 1e-10)
        fprintf("%3d %10.3e \n", i, norm(r));
        break;
    end % End if 
    
%     fprintf("%3d   %10.3e %10.3e \n", i, alpha, norm(r));
    
end % End for

fprintf("%3d %10.3e \n", i, norm(r));

end % End function