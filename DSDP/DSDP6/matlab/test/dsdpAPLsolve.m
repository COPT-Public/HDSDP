function [yut, fub, iter] = dsdpAPLsolve(y0, A, C, iter, theta, lbd)
% Implement the accelerated prox-level algorithm for LMI problem

if nargin < 6
    lbd = 0.0;
end % End if

if nargin < 5
    theta = 0.5;
end % End if 

[fu0, ~] = eigrad(A, C, y0);

fub = fu0;
yt = y0;
yut = y0;

for i = 1:iter
    
    alpha = 2 / (i + 1);
    ylt = (1 - alpha) * yut + alpha * yt;
    
    [fylt, grad] = eigrad(A, C, ylt);
    
    if fylt + grad' * (yt - ylt) >= -lbd
        yt = yt - (fylt + lbd) / norm(grad)^2 * grad;
    end % End if
    
    ytmp = alpha * yt + (1 - alpha) * yut;
    
    [fytmp, ~] = eigrad(A, C, ytmp);
    
    % fprintf("%f \n", fytmp);
    
    if fytmp < fub
        yut = ytmp;
        fub = fytmp;
    end % End if
    
    if (fub + lbd) < theta * (fu0 + lbd)
        break;
    end % End if
    
end % End for

iter = i;

end % End function