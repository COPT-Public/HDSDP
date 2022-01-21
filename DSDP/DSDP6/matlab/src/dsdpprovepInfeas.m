function [y, S, reason] = dsdpprovepInfeas(A, b, C, y, S, kappa, mu, reason, dsdpParam)
% Prove primal infeasibility of the problem
maxiter = dsdpParam{24};
stepstrategy = dsdpParam{15};
ndash        = dsdpParam{20};
Ry = - C; 

m = length(b);
[n, ~] = size(C);

% Ry = - speye(n, n);
% S = speye(n, n);
% y = zeros(n, 1);

fprintf("%4s  %10s  %10s  %10s \n", "Iter", "dInf", "dObj", "step");
showdash(ndash);
reason = "DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE";

% for i = 1:maxiter
%     
%     [M, ~, asinv, ~, ~, ~, ~, ~, ~, ~] = dsdpgetSchur(A, S);
%     
%     if mu < 1e-05
%         M = M + eye(m, m) * 1e-08;
%     end % End if 
%     
%     dy1dy2 = M \ [b, asinv];
%     dy1 = dy1dy2(:, 1) / tau;
%     dy2 = dy1dy2(:, 2);
%     
%     dymu = dy1 / mu - dy2;
%     newmu = dsdpselectMu(A, S, muprimal, dymuprimal, dy1, backwardnewton, ismufeas);
%     delta = sqrt(dymu' * M * dymu);
%     
%     if delta < 0.1
%         mu = mu * 0.1;
%     end % End if 
%     
%     
%     
% end % End for

for i = 1:maxiter
    
    dObj = b' * y;
    [M, u, asinv, ~, ~, ~, csinv, ~, ~, ~] = dsdpgetSchur(A, S, C);
    
    d3 = M \ asinv;
    
    dy = -2 * d3 - y;
    dkappa = (b - mu * u)' * (y + dy) - 2 * mu * csinv + dObj - kappa;
    dS = Ry - dsdpgetATy(A, dy);
    
    try
        step = dsdpgetStepsize(S, dS, kappa, dkappa, 1.0, 1.0, stepstrategy, 0.95);
    catch
        return;
    end % End try
    
    y = y + step * dy;
    kappa = kappa + step * dkappa;
    Ry = Ry * (1 - step);
    S = - Ry - dsdpgetATy(A, y);
    
    y = y / kappa;
    S = S / kappa;
    Ry = Ry / kappa;
    kappa = 1;
    
    nrmRd = norm(Ry, 'fro');
    
    fprintf("%4d  %10.3e  %10.3e  %10.3e \n", i, nrmRd, dObj, step);
    if step < 1e-05
        break;
    end % End if
    
    
    if nrmRd / dObj < 1e-10 && dObj > 0
        showdash(ndash);
        fprintf("Dual ray found \n");
        reason = "DSDP_PRIMAL_INFEASIBLE_DUAL_UNBOUNDED";
        showdash(ndash);
        break;
    end % End if 
    
    % mu = mu * 0.1;
    
end % End if

end % End function