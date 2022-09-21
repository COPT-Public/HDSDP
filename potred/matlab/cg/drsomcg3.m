function [x, e] = drsomcg3(A, b, x0, tol)
% Solve linear system using dimension-reduced second order method
n = size(A, 1);
if nargin >= 3
    x_prev = x0;
else
    tol = 1e-10;
    x_prev = zeros(n, 1);
end % End if

r_prev = A * x_prev - b;
alpha = r_prev' * r_prev; alpha = alpha / (r_prev' * A * r_prev);
x_pres = x_prev - alpha * r_prev;

G = eye(3);
warning off;
% Minimize the quadratic form 0.5 * x' * A * x - b' * x
for i = 1:100
    
    gk = A * x_pres - b;
    dk = x_pres - x_prev;    
    mk = A * dk;
    nrmgk = norm(gk);
    
    if nrmgk < tol
        break;
    end % End if 
    
    Agk = A * gk; Adk = A * dk; Amk = A * mk;
    
    gkAgk = gk' * Agk; 
    gkAdk = dk' * Agk; dkAdk = dk' * Adk; 
    gkAmk = mk' * Agk; dkAmk = dk' * Amk; mkAmk = mk' * Amk;
    
    gkgk = norm(gk)^2;
    gkdk = gk' * dk;
    gkmk = gk' * mk;
    
    Q = [ gkAgk, -gkAdk, -gkAmk;
         -gkAdk,  dkAdk,  dkAmk;
         -gkAmk,  dkAmk,  mkAmk];
     
%     G = [ gkgk, -gkdk, -gkmk;
%          -gkdk,  norm(dk)^2,  dkmk;
%          -gkmk,  dkmk,  norm(mk)^2];
     
    r = 0.0;
    Q = Q + r * G;
    
    c = [-gkgk; gkdk; gkmk];
     
    alpha = - Q \ c;
    
%     assert(norm(Q * alpha - c) < 1e-10); 
    
    x_prev = x_pres;
    x_pres = x_pres - alpha(1) * gk + alpha(2) * dk + alpha(3) * mk;
    
end % End for

x = x_pres;
e = gk;

fprintf("%3d %10.3e %10.3e\n", i, fx(A, b, x_pres), nrmgk);

end % End function