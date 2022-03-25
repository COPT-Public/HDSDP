function [S, y, kappa, tau, mu, pObj, reason, iter] = dsdpDualInfeasPhaseRlx(A, b, C, y, S, Rd, mu, kappa, tau, dsdpParam)
% Implementation of dual scaling algorithm phase 1: dual infeasibility
% elimination

[n, ~] = size(S);

maxiter      = dsdpParam{1};
tol          = dsdpParam{2};
nrmC         = norm(C, 'fro');
alphaphase1  = dsdpParam{14};
sigma        = dsdpParam{3};
initstrategy = dsdpParam{7};
stepstrategy = dsdpParam{15};
ncorrp1      = dsdpParam{9};
corralpha    = dsdpParam{16};
candmu       = dsdpParam{19};
ndash        = dsdpParam{20};
corrdelta    = dsdpParam{22};
m = length(y);
eyem = eye(m) * 1e-10;
pObj = inf;
muprimal = mu;
augment = false;
nrmrtk = 0;
step = 0;
delta = inf;

for i = 1:maxiter
    
    if initstrategy == "IdS"
        nrmRd = norm(Rd, 'fro') / (tau * (1 + nrmC));
    else
        nrmRd = sqrt(n) * abs(full(Rd(1, 1))) / (tau * (1 + nrmC));
    end % End if
    
    if nrmRd < tol / 10
        if ~ isinf(pObj)
            reason = "DSDP_PRIMAL_DUAL_FEASIBLE";
        else
            reason = "DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE";
        end % End if 
        break;
    end % End if  
    
    if tau < (0.001 * kappa) && mu < tol
        if ~ isinf(pObj)
            reason = "DSDP_PRIMAL_FEASIBLE_DUAL_INFEASIBLE";
        else
            reason = "DSDP_DUAL_INFEASIBLE";
        end % End if 
        break;
    end % End if
    
    dObj = b' * y;
    nrmrtk = abs(tau * kappa - mu);
    
    [M, u, asinv, ~, ~, ~, csinv, csinvcsinv, asinvrysinv, csinvrysinv] = ...
        dsdpgetSchur(A, S, C, Rd, initstrategy);
    
    d2_12_3_4 = M \ [b, u, asinv, asinvrysinv];
    
    d2  = d2_12_3_4(:, 1);
    d12 = d2_12_3_4(:, 2);
    d3  = d2_12_3_4(:, 3);
    d4  = d2_12_3_4(:, 4);
       
    % Check primal feasibility
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for newmu = candmu * muprimal
        [delta, newpObj] = dsdpgetProxMeasure(d12, d2, d3, u, csinvcsinv, csinv, tau, newmu, dObj, M, A, C, y, Rd);
        if ~ isinf(newpObj)
            pObj = newpObj / tau;
            if newmu > tol^2
                muprimal = newmu;
            end % End if 
            mu = max(muprimal, mu * sigma);
            break;
        end % End if
    end % End for
    
    b1 = b - mu * u;
    b2 = d2 * tau / mu - d3 + d4;
    d11 = d2 / mu;
    d1  = d11 + d12;
       
    % Take Newton step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    taudenom = b1' * d1 + mu * csinvcsinv + kappa / tau;
    if abs(taudenom) < tol^3
        dtau = 0;
    else
        dtau = - b' * y + mu * (1 / tau + csinv - csinvrysinv) - b1' * b2;
        dtau = dtau / taudenom;
    end % End if
    
    dy     = d1 * dtau + b2;
    dS     = Rd + C * dtau - dsdpgetATy(A, dy);
    dkappa = - kappa + (mu / tau) - (kappa / tau) * dtau;
    
    step = dsdpgetStepsize(S, dS, kappa, dkappa, tau, dtau, stepstrategy, alphaphase1);
    
    if step < 1e-03
        reason = "DSDP_SMALL_STEP";
        break;
    end % End if 
    
    y = y + step * dy;
    % S = S + step * dS;
    kappa = kappa + step * dkappa;
    tau = tau + step * dtau;
    Rd = Rd * (1 - step);
    
    if step < 0.5 && augment
        augment = false;
        Rd = Rd - speye(n, n);
    end
    
    S = - dsdpgetATy(A, y) + C * tau - Rd;
    
    % assert(pObj >= dObj / tau);
    % assert(min(eig(S)) > 0);
    
    % Corrector 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [y, S] = dinfeaspotrdc(A, b * tau, C * tau, y, S, Rd, M, d2 * tau, mu, ncorrp1);
    % [y, S, kappa, tau] = dsdpCorrector(A, b, C, y, S, Rd, kappa, tau, mu, M, d2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Logging
    fprintf("%3d  %10.2e  %10.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e\n",...
        i, pObj, b' * y / tau, nrmRd, nrmrtk, muprimal, mu, step, delta);
    
    % TODO: decrese mu
    if delta < 0.1
        mu = mu * 0.1;
    end % End if 
    
    if i > 30
        sigma = 0.1;
        alphaphase1 = 0.2;
    end % End if 
    
    nrm = norm(y, 'inf');
    if nrm > 1e+08
        y   = y / nrm;
        S   = S / nrm;
        tau = tau / nrm;
        Rd  = Rd / nrm;
    end % End if 
    
end % End for

% [y, S] = dinfeaspotrdc(A, b, C * tau, y, S, sparse(n, n), M, d2 * tau, mu, 12);

iter = i;

% Logging
fprintf("%3d  %10.2e  %10.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e \n",...
        i, pObj, b' * y / tau, nrmRd, nrmrtk, muprimal, mu, step, delta);
showdash(ndash);
if reason == "DSDP_DUAL_INFEASIBLE"
    fprintf("Phase 1 certificates dual infeasibility \n");
elseif reason == "DSDP_SMALL_STEP"
    fprintf("Phase 1 ends due to small stepsize \n");
elseif reason == "DSDP_DUAL_FEASIBLE"
    fprintf("Phase 1 finds a dual feasible solution \n");
end % End if
    
end % End function