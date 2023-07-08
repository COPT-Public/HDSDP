function [x, y, s] = dualpot(A, b, c, verbose)
% Solve LP using two-phase dual potential reduction
% This part of the code solves
% min c' * x
%      A * x == b
%          x >= 0
if nargin == 3
    verbose = 0;
end % End if

warning off;
[m, n] = size(A);
x = ones(n, 1);
y = zeros(m, 1);

nrmb = norm(b);
nrmc = norm(c);

if nrmb > 1e+05
    nrmb = 1e+04;
end % End if 

if nrmc > 1e+05
    nrmc = 1e+04;
end % End if 

[row_id, col_id, v_elem] = find(A);

b = b / nrmb;
c = c / nrmc;

tau = 1;
Rd = - norm(c) * 10;
s = c * tau - Rd * tau;
mu = 1.0;

if verbose
    fprintf("%4s  %5s  %8s  %5s\n", "Iter", "dRes", "dObj", "mu");
end % End if 

for i = 1:100
    
    % Set up the Schur complement
    ASinv = sparse(row_id, col_id, v_elem ./ s(col_id));
    M = ASinv * ASinv';
    
    sinv = 1 ./ s;
    ASinv = A * sinv;
    ASinvRdSinv = Rd * (A * (sinv.^2));
    ASinvCSinv = A * ((sinv.^2) .* c);
    CSinvRdSinv = Rd * (c' * (sinv.^2));
    CSinv = c' * sinv;
    CSinvCSinv = norm(c .* sinv)^2;
    
    % Set up the KKT system
    % In the HSD model, we first solve
    % M * d1 = b 
    % M * d2 = ASinv
    % M * d3 = ASinvRdSinv
    % M * d4 = ASinvCSinv   => Only solved in HSD
    try
        Mdcp = decomposition(M, 'chol');
    catch
        break;
    end % End try
    
    d1 = Mdcp \ b;
    d2 = Mdcp \ ASinv;
    d3 = Mdcp \ ASinvRdSinv;
    d4 = Mdcp \ ASinvCSinv;
    
    dd1 = b - mu * ASinvCSinv;
    
    dObj = b' * y;
    dtauUp = -dObj + mu / tau + mu * (CSinv - CSinvRdSinv) - dd1' * (d1 * tau / mu - d2 + d3);
    dtauLow = mu * CSinvCSinv + mu / tau^2 + dd1' * (d1 / mu + d4);
    
    if abs(dtauLow) < 1e-12
        dtau = 0;
    else
        dtau = dtauUp / dtauLow;
    end % End if
    
    dy = d1 * (tau + dtau) / mu + d4 * dtau - d2 + d3;
    ds = Rd - A' * dy + c * dtau;
    
    alphas = 1 / min(ds./s);
    alphatau = tau / dtau;
    
    if alphatau > 0
        alphatau = 100;
    end % End if 
    
    if alphas > 0
        alphas = 100;
    end % End if
    
    if alphas > 0 && alphatau > 0
        alpha = 1.0;
    else
        alpha = min([1.0, 0.9 * min(abs([alphas, alphatau]))]);
    end % End if
    
    y = y + alpha * dy;
    tau = tau + alpha * dtau;
    Rd = Rd * (1 - alpha);
    s = - Rd - A' * y + c * tau;
    
    if verbose
        fprintf("%3d  %5.2e  %8.5e  %5.2e\n", i, abs(Rd / tau), b' * y / tau * nrmb * nrmc, mu);
    end % End if 
    
    if Rd ~= 0
        mu = min(max(mu * 0.7, abs(Rd / tau) * 0.1), mu);
    else
        mu = mu * 0.1;
    end % End if
    
%     if min(s) <= 1e-14
%         scal_factor = 10;
%         y = y * scal_factor;
%         s = s * scal_factor;
%         tau = tau * scal_factor;
%         Rd = Rd * scal_factor;
%         mu = mu * scal_factor^2;
%     end % End if
    
    if mu < 1e-08 && abs(Rd / tau) < 1e-08
        break;
    end % End if
    
end % End for

y = y / tau * nrmc;
s = s / tau * nrmc;
x = mu ./ s;

end % End function