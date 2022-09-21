function [x, e] = drsomcg(A, b, x0, tol)
% Solve linear system using dimension-reduced second order method
n = size(A, 1);
if nargin >= 3
    x_prev = x0;
else
    tol = 1e-10;
    x_prev = zeros(n, 1);
end % End if


r = A * x_prev - b;
alpha = r' * r; alpha = alpha / (r' * A * r);
x_pres = x_prev - alpha * r;

% trA = gallery('tridiag', n, 1, 1, 1);
% trA = A .* trA;
xstar = A \ b;

% Minimize the quadratic form 0.5 * x' * A * x - b' * x
for i = 1:100
    
    gk = A * x_pres - b;
    
    % The momentum direction
    if i <= 10
        dk = 1e-07 * randn(n, 1) + (xstar - x_pres);
    else
        dk = (x_pres - x_prev);
    end % End if
    
%     dk = (gk ./ diagA) - x_pres;
    
    rnrm = norm(gk);
    
    if rnrm < tol
        break;
    end % End if 
    
    Agk = A * gk; Adk = A * dk;
    gkAgk = gk' * Agk; dkAgk = dk' * Agk; dkAdk = dk' * Adk;
    
    gkdk  = gk' * dk;
    
    Q = [ gkAgk,  -dkAgk;
         -dkAgk,  dkAdk];
    c = [-rnrm^2; gkdk];
    
    if rnrm < 1e-12
        s = norm(Q, 'fro');
        Q = Q / s;
        c = c / s;
    end % End if
          
    
%     alpha = - Q \ c;
     
    alpha = - solveQ(Q, c);
    
%     if i >= 70
%         ;
%     end % End if
        
    x_prev = x_pres;
    x_pres = x_pres - alpha(1) * gk + alpha(2) * dk;
    
    fprintf("%3d %10.3e %10.3e\n", i, fx(A, b, x_pres), rnrm);
    
end % End for

x = x_pres;
e = gk;

fprintf("%3d %10.3e %10.3e\n", i, fx(A, b, x_pres), rnrm);

end % End function