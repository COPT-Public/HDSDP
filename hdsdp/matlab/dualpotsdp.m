function [X, y, S] = dualpotsdp(A, b, C, verbose)
% Solve LP using two-phase dual potential reduction
% This part of the code solves
% min c' * x
%      A * x == b
%          x >= 0
if nargin == 3
    verbose = 0;
end % End if

warning off;
[m, ~] = size(A);
y = zeros(m, 1);

n = size(C, 1);

tau = 1;
Rd = - 100 * norm(C, 'fro');
S = C * tau - eye(n) * Rd * tau;
mu = 1e+08;

if verbose
    fprintf("%4s  %5s  %8s  %5s\n", "Iter", "dRes", "dObj", "mu");
end % End if 

for i = 1:100
    
    % Set up the Schur complement
    [M, ASinv, ASinvRdSinv, ASinvCSinv, CSinv, CSinvCSinv, CSinvRdSinv] =...
        hdsdp_kktbuild(A, C, S, Rd);
    
    % Set up the KKT system
    % In the HSD model, we first solve
    % M * d1 = b 
    % M * d2 = ASinv
    % M * d3 = ASinvRdSinv
    % M * d4 = ASinvCSinv   => Only solved in HSD
    try
        M = M + eye(m) * 1e-10;
        Mdcp = decomposition(M);
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
    
    if abs(dtauLow) < 1e-14
        dtau = 0;
    else
        dtau = dtauUp / dtauLow;
    end % End if
    
    dy = d1 * (tau + dtau) / mu + d4 * dtau - d2 + d3;
    dS = Rd * eye(n) - hdsdp_aty(A, dy) + C * dtau;
    
    alphas = hdsdp_ratiotest(S, dS);
    alphatau = tau / dtau;
    
    if alphatau > 0
        alphatau = 100;
    end % End if 
    
    alpha = min([1.0, 0.7 * min(abs([alphas, alphatau]))]);
    
    y = y + alpha * dy;
    tau = tau + alpha * dtau;
    Rd = Rd * (1 - alpha);
    S = - Rd * eye(n) - hdsdp_aty(A, y) + C * tau;
    
    try 
        chol(S);
    catch
        keyboard;
    end % End try
    
    if verbose
        fprintf("%3d  %5.2e  %+10.8e  %5.2e tau: %5.2e alpha: %4.3f\n", i, abs(Rd / tau) * sqrt(n), b' * y / tau, mu, tau, alpha);
    end % End if 
   
    if alpha > 0.6 && tau > 1
        mu = min(max(mu * 0.1, abs(Rd / tau) * 0.1), mu);
%     elseif alpha < 0.5 || tau < 0.1
%         mu = 0.9 * mu;
    else
        mu = min(max(mu * 0.5, abs(Rd / tau) * 0.1), mu);
    end % End if
    
    if mu < 1e-08 && abs(Rd / tau) < 1e-08
        break;
    end % End if
    
    if tau < 1e-08
        fprintf("Suspected dual infeasible \n");
        break;
    end % End if
    
end % End for

y = y / tau;
S = S / tau;
X = mu * inv(S);

end % End function