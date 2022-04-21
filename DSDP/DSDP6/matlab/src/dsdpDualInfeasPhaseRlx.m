function [S, y, kappa, tau, mu, pObj, reason, iter] = dsdpDualInfeasPhaseRlx(A, b, C, y, S, Rd, mu, kappa, tau, dsdpParam)
% Implementation of dual scaling algorithm phase 1: dual infeasibility
% elimination

[n, ~] = size(S);
m = length(b);
maxiter      = dsdpParam{1};
tol          = dsdpParam{2};
nrmC         = norm(C, 'fro');
alphaphase1  = dsdpParam{14};
sigma        = dsdpParam{3};
initstrategy = dsdpParam{7};
stepstrategy = dsdpParam{15};
ncorrp1      = dsdpParam{9};
candmu       = dsdpParam{19};
ndash        = dsdpParam{20};

pObj = inf;
muprimal = mu;
step = 0;
delta = inf;
bound = 1e+07;
ub = bound;
lb = -ub;
sl = y - lb * tau;
su = ub * tau - y;

pweight      = 0.0;% dsdpParam{29};
prelax       = true;

if prelax
    np = 2 * m + n;
else
    np = n;
end % End if 

corrrhs = zeros(m, 1);

for i = 1:maxiter
    
    nrmRd = sqrt(n) * abs(full(Rd(1, 1))) / (tau * (1 + nrmC));
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
    
    [M, u, asinv, L, ~, ~, csinv, csinvcsinv, asinvrysinv, csinvrysinv, rysinv] = ...
        dsdpgetSchur(A, S, C, Rd, initstrategy);
    
    % Primal relaxation
    if prelax
        M = M + diag(sl.^-2 + su.^-2);
        u = u + bound * (sl.^-2 + su.^-2);
        asinv = asinv - sl.^-1 + su.^-1;
        csinv = csinv + bound * (sum(su.^-1) - sum(sl.^-1));
        csinvcsinv = csinvcsinv + (bound * bound) * sum(sl.^-2 + su.^-2);
    end % End if 
    d2_12_3_4 = M \ [b, u, asinv, asinvrysinv];
    
    d2  = d2_12_3_4(:, 1);
    d12 = d2_12_3_4(:, 2);
    d3  = d2_12_3_4(:, 3);
    d4  = d2_12_3_4(:, 4);
    
    % Check primal feasibility
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for newmu = mu
        [delta, newpObj] = dsdpgetProxMeasure(d2, d3, tau, newmu, dObj, M,...
            A, C, y, Rd, asinvrysinv, rysinv, asinv, S, b, np, bound);
        delta = sqrt(delta / tau^2);
        if ~ isinf(newpObj)
            if ~isinf(pObj)
                pObj = min(pObj, newpObj);
                muprimal = (pObj - dObj - Rd(1, 1) / tau * 1e+06) / (np * 3);
                mu = min(mu, muprimal);
            else
                pObj = newpObj / tau;
                muprimal = (pObj - dObj - Rd(1, 1) / tau * 1e+06) / (np * 3);
                mu = muprimal;
            end % End if
        end % End if
        
        break;
    end % End for

    if delta < 0.1
        mu = mu * 0.1;
    end % End if
    
    % b1  = b * pweight - mu * u;
    b2  = d2 * (pweight * tau / mu) - d3 + d4;
    d11 = d2 / mu;
    d1  = pweight * d11 + d12;
    
    % Adaptive drate
    if (true)
        dSd3 = dsdpgetATy(A, d3);
        alphac = dsdpgetalpha(S, dSd3);
        
        if alphac < 0
            alphac = 1;
        else
            alphac = min(alphac * 0.95, 1.0);
        end % End if 
        
        barrierold = 2 * sum(log(diag(L)));
        
        while alphac > 1e-04
            Shat = S + alphac * dSd3;
            Lhat = lchol(Shat);
            barriernew = 2 * sum(log(diag(Lhat))); 
            if barriernew > barrierold || delta < 10.0 || true
                break;
            else
                alphac = alphac * 0.6;
            end % End if 
        end % End while

        gammafeas = dsdpgetalpha(Shat, Rd - dsdpgetATy(A, alphac * d4));
        drate = min(gammafeas * 0.95, 1.0);
        b2 = -d3 + drate * d4;
    end % End if
    
    % Compute Newton step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     taudenom = b1' * d1 + mu * csinvcsinv + kappa / tau;
%     if abs(taudenom) < 1e-15
%         dtau = 0;
%     else
%         dtau = - pweight * b' * y + mu * (1 / tau + csinv - csinvrysinv) - b1' * b2;
%         dtau = dtau / taudenom;
%     end % End if
    
    dtau = 0.0;
    dy  = d1 * dtau + b2;
    dS  = drate * Rd + C * dtau - dsdpgetATy(A, dy);
    
    if prelax
        dsu = -dy + ub * dtau;
        dsl =  dy - lb * dtau;
    end % End if 
    
    % Compute stepsize
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    step = dsdpgetStepsize(S, dS, kappa, 1.0, tau, dtau, stepstrategy, 1.0);
    
    if prelax
        tmp = - 1 / min(dsl ./ sl);
        if tmp > 0
            step = min(step, tmp);
        end % End if
        tmp = - 1 / min(dsu ./ su);
        if tmp > 0
            step = min(step, tmp);
        end % End if
    end % End if
    
    refstep = step;
    step = min(step * 0.95, alphac);
    
    if step < 1e-03 && (false)
        reason = "DSDP_SMALL_STEP"; 
        break;
    end % End if
    
    % Take Newton step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = y + step * dy;
    tau = tau + step * dtau;
    kappa = mu / tau;
    Rd = Rd * (1 - drate * step);
    S = - dsdpgetATy(A, y) + C * tau - Rd;
    
    % Dual infeasibility heuristic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if refstep > 1
%         if drate > 0.5
%             ncorrp1 = floor(ncorrp1 / 2);
%         end % End if 
%     end % End if 
%     
%     if refstep < 0.1
%         ncorrp1 = min(12, 2 * max(ncorrp1, 2));
%     elseif refstep < 0.5
%         ncorrp1 = min(8, 2 * max(ncorrp1, 1));
%     end % End if 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    assert(min(y) > -1e+07);
    assert(max(y) < 1e+07);
    
    % Corrector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [y, S, Rd] = dinfeaspotrdc(A, b * tau, C * tau, y, S, Rd, M, d2 * tau, mu, 12, bound, 0.9);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    assert(min(y) > -1e+07);
    assert(max(y) < 1e+07);
    
    if prelax
        sl = y - lb * tau;
        su = ub * tau - y;
    end % End if 
    
    % Logging
    dObj = b' * y;
    fprintf("%3d  %10.2e  %10.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e\n",...
        i, pObj, dObj / tau, nrmRd, drate, mu, step, delta);
    
end % End for

% Corrector at the end
[y, S] = dinfeaspotrdc(A, b, C * tau, y, S, sparse(n, n), M, d2 * tau, mu, 4, bound, 0.0);

iter = i;

% Logging
fprintf("%3d  %10.2e  %10.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e \n",...
    i, pObj, b' * y / tau, nrmRd, kappa / tau, muprimal, step, delta);
showdash(ndash);
if reason == "DSDP_DUAL_INFEASIBLE"
    fprintf("Phase 1 certificates dual infeasibility \n");
elseif reason == "DSDP_SMALL_STEP"
    fprintf("Phase 1 ends due to small stepsize \n");
elseif reason == "DSDP_DUAL_FEASIBLE"
    fprintf("Phase 1 finds a dual feasible solution \n");
end % End if

end % End function

%     % Complicated corrector
%     if false
%         invS = S \ speye(n);
%         SinvdSSinvdS = invS * dS * invS * dS * invS;
%         for k = 1:m
%             corrrhs(k) = trace(A{k} * SinvdSSinvdS);
%         end % End for
%         dycorr = - (M \ corrrhs);
%         dy = dy + dycorr;
%         dS = drate * Rd - dsdpgetATy(A, dy);
%     end % End if