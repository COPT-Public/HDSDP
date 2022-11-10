function [x, y, s, kappa, tau] = hsdipm(A, b, c)
% Implement a primal-dual interior point method using HSD embedding

[m, n] = size(A);
x = ones(n, 1);
y = zeros(m, 1);
s = ones(n, 1);
tau = 1;
kappa = 1;

sigma = 0.7;

fprintf("%8s %8s %8s %8s %8s \n", "pobj", "dobj", "pinf", "dinf", "mu");

for i = 1:30
    
    mu = (x' * s + kappa * tau) / (n + 1);
    
    % Set up residuals
    rp = A * x - b * tau;
    rd = -A' * y - s + c * tau;
    pobj = c' * x;
    dobj = b' * y;
    
    fprintf("%8.1e %8.1e %8.1e %8.1e %8.1e \n", pobj / tau, dobj / tau, norm(rp) / 1, norm(rd) / 1, mu);
    
    XSe = sqrt(x .* s);
    D = sqrt(s) ./ sqrt(x);
    Dinv = sqrt(x) ./ sqrt(s);
    
    DinvXinvrmu1 = XSe - (mu * sigma) ./ XSe;
    
    ADinv = A * diag(Dinv);
    
    M = [speye(n), ADinv';
         ADinv, sparse(m, m)];
    
    Dinvc = c ./ D;
    rhs1 = [-Dinvc; b];
    rhs2 = [Dinvc; b];
    
    d1 = M \ rhs1;
    d2 = M \ [rd ./ D + DinvXinvrmu1; rp];
    dtau = rhs2' * d2 + dobj - pobj - sigma * mu / tau;
    dtau = - dtau / (kappa / tau - rhs2' * d1);
    
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
    
end % End for

end % End function