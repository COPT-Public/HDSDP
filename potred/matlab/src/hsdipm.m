function [x, y, s, kappa, tau] = hsdipm(A, b, c)
% Implement a primal-dual interior point method using HSD embedding

[m, n] = size(A);
x = ones(n, 1);
y = zeros(m, 1);
s = ones(n, 1);
tau = 1;
kappa = 1;

sigma = 0.7;

% fprintf("%4s %8s %8s %8s %8s %8s \n", "iter", "pobj", "dobj", "pinf", "dinf", "mu");

for i = 1:100
    
    mu = (x' * s + kappa * tau) / (n + 1);
    
    % Set up residuals
    rp = A * x - b * tau;
    rd = -A' * y - s + c * tau;
    pobj = c' * x;
    dobj = b' * y;
    
    pinf = norm(rp) / tau;
    dinf = norm(rd) / tau;
    
    if max([pinf, dinf, mu]) < 1e-06
        break;
    end % End if
    
%     fprintf("%4d %8.1e %8.1e %8.1e %8.1e %8.1e \n", i, pobj / tau, dobj / tau, pinf, dinf, mu);
    
    XSe = sqrt(x .* s); % n 
    D = sqrt(s) ./ sqrt(x); % n
    Dinv = sqrt(x) ./ sqrt(s);
    
    ADinv = A * diag(Dinv);
    
    M = [speye(n), ADinv';
         ADinv, sparse(m, m)];
    
    Dinvc = c ./ D; 
    DinvXinvrmu1 = XSe - (mu * sigma) ./ XSe; 
    rhs1 = [-Dinvc; b]; % m + n
    rhs2 = [rd ./ D + DinvXinvrmu1; rp];
    aux = [Dinvc; b]; % m + n
    
    d1 = M \ rhs1; 
    d2 = M \ rhs2;
    
    dtau = aux' * d2 + dobj - pobj - sigma * mu / tau; 
    dtau = - dtau / (kappa / tau - aux' * d1);
    
    dxdy = d1 * dtau - d2;
    dx = dxdy(1:n) ./ D;
    dy = - dxdy(n+1:end);
    dkappa = - kappa * dtau / tau - kappa + sigma * mu / tau;
    ds = -D.* dxdy(1:n) - s + (mu * sigma) ./ x;
    
    alpha = 0.995 / abs(min([dx./x; ds./s; dkappa / kappa; dtau./tau]));
    
    x = x + alpha * dx;
    y = y + alpha * dy;
    s = s + alpha * ds;
    kappa = kappa + alpha * dkappa;
    tau = tau + alpha * dtau;
    
    simp = tau; % sum(x) + sum(s) + kappa + tau;
    x = x / simp;
    y = y / simp;
    s = s / simp;
    kappa = kappa / simp;
    tau = tau / simp;
    
end % End for

end % End function