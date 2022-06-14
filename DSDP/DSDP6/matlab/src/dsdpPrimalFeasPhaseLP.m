function [y, S, pObj, muprimal, xmaker, reason] = dsdpPrimalFeasPhaseLP(A, b, C, y, S, muHSD, pObj, dObj, tau, kappa, ispdfeas, dsdpParam, bound)
% Implementing primal feasibility certificate for DSDP
% y and S are input such that ATy + S - C * tau = 0
% As if using original DSDP to solve
% maximize tau * b' * y s.t. ATy + S - C * tau = 0
% Assume that we now have bound constraint l <= y <= u
[n, ~] = size(S);
m = length(y);
rhouser      = dsdpParam{6};
maxiter      = dsdpParam{1};
maxpfeasiter = dsdpParam{17};
ncorr        = dsdpParam{27};
tol          = dsdpParam{2};
ndash        = dsdpParam{20};
tau          = 1;
normalizer   = dsdpParam{23};
pinfeasbound = dsdpParam{26};
dObj = b' * y;
isfirst      = true;

pinfeas = inf;

if ispdfeas
    maxpfeasiter = dsdpParam{1};
end % End if

xmaker = cell(3, 1);
smallstep = zeros(3, 1);

reason = "DSDP_MAXITER";

u = bound;
l = -u;
sl = y - l;
su = u - y;

muprimal = min((pObj - dObj) / rhouser, muHSD);
bettermu = false;

for i = 1:maxiter
    
    if i >= maxpfeasiter + 1 && reason == "DSDP_MAXITER"
        reason = "DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE";
        break;
    end % End if
    
    ismufeas = false;
    muHSD = min(muHSD, muprimal);
    
    % Schur of SDP
    [M, ~, asinv, ~, ~, ~, ~, ~, ~, ~] = dsdpgetSchur(A, S);
    
    % Schur of LP
    slinv2 = sl.^-2;
    suinv2 = su.^-2;
    asinv = asinv - sl.^-1 + su.^-1;
    M = M + diag(slinv2) + diag(suinv2);
    
    diagM = max(diag(M)); 
    Mhat = M;
    
    if muprimal < 1e-05
        Mhat = M + eye(m, m) * 1e-13;
    end % End if
    
    if isfirst
        Mhat = Mhat + eye(m, m) * min(1e-12, diagM * 1e-08);
    end % End if
        
    % dy1 = dsdpConjGrad(Mhat, b, diagM, zeros(m, 1)) / tau;
    % dy2 = dsdpConjGrad(Mhat, asinv, diagM, zeros(m, 1));
    dy1dy2 = Mhat \ [b, asinv];
    dy1 = dy1dy2(:, 1);
    dy2 = dy1dy2(:, 2);
    
    dymuprimal = dy1 / muprimal - dy2;
    backwardnewton = C - dsdpgetATy(A, y - dymuprimal);
    backwardnewtonub = u - (y - dymuprimal);
    backwardnewtonlb = -l + (y - dymuprimal);
    
    % Proximity
    delta = sqrt(dymuprimal' * (b / muprimal - asinv));
    if (true)
        csinv = trace(S \ C); %#ok
        [~, pObjtmp, ~] = dsdpgetmualpha(asinv, b, csinv, muprimal, bound, 1e-06);
        
        if pObjtmp < pObj
            pObj = pObjtmp;
        end
    end % End if
    
    % Check feasibility of backward Newton step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if dsdpIspsd(backwardnewton) && min(backwardnewtonlb) >= 0 && min(backwardnewtonub) >= 0
        ismufeas = true;
        reason = "DSDP_PRIMAL_DUAL_FEASIBLE";
        diff = muprimal * (dymuprimal' * asinv + (n + 2 * m));
        
        if (false)
            Xtmp = muprimal * (S \ (S \ (backwardnewton))'); %#ok
            xl = muprimal * (sl.^-2).*backwardnewtonlb;
            xu = muprimal * (su.^-2).*backwardnewtonub;
            ax = zeros(m, 1);
            for q = 1:m
                ax(q) = trace(A{q} * Xtmp);
            end % End for
            resi = ax - xl + xu - b;
        end % End if
        
        pObj = min(dObj + diff, pObj);
        muub = (pObj - dObj) / (n + 2 * m);
        mulb = muub / rhouser;
        newub = "*";
        xmaker{1} = y;
        xmaker{2} = dymuprimal;
        xmaker{3} = muprimal;
        % Get primal infeasibility
        pinfeasl = muprimal * norm(sl.^-1 + slinv2 .* dymuprimal, 'inf'); 
        pinfeasu = muprimal * norm(su.^-1 - suinv2 .* dymuprimal, 'inf'); 
        pinfeas = max(pinfeasl, pinfeasu);
        % fprintf("pinfeas: %10.3e \n", pinfeas);
    else
        % fprintf("pinf: %e pObj: %e \n", norm(infeastmp), mutmp * trace(S \ C));
        muub = (pObj - dObj) / (n + 2 * m);
        mulb = muub / rhouser;
        newub = "";
    end % End if
    
    % Select new mu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newmu = dsdpselectMuRlx(A, sl, su, S, muprimal, dymuprimal, dy1, backwardnewton, backwardnewtonlb, backwardnewtonub, ismufeas, (pObj - dObj) / (n + 2 * m));
    muprimal = min(newmu, muub);
    muprimal = max(muprimal, mulb);
    
    if delta < 0.1 && ismufeas
        muprimal = muprimal * 0.1;
    end % End if
    
    % Potential reduction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dymuprimal = dy1 / muprimal - dy2;
    dS = - dsdpgetATy(A, dymuprimal);
    dsu = - dymuprimal;
    dsl = dymuprimal;
    
    [y, S, sl, su, step] = dsdptakedualStepRlx(A, b, C, l, u, pObj, y, dymuprimal, S, dS, sl, dsl, su, dsu, (pObj - dObj) / muprimal);
    % [y, S, step] = dsdptakedualStep(A, b, C, pObj, y, dymuprimal, S, dS, rhouser * (n + sqrt(n)));
    assert(dsdpIspsd(S));
    
    if step < 1e-05 && ~ismufeas
        if sum(smallstep) == 1
            reason = "DSDP_SMALL_STEP";
            break;
        else
            [~, idx] = min(smallstep);
            smallstep(idx) = 1;
        end % End if
    else
        smallstep = zeros(3, 1);
    end % End if
    
    % Corrector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ncorr > 0
        if delta < 0.1
            ncorr = 0;
        end % End if
        
        if step < 0.1 && muprimal < 1e-05
            ncorr = 0;
        end % End if
        
        if muprimal < 1e-06
            ncorr = 0;
        end % End if
    end % End if
    
    [y, S, muprimal] = dsdpdualCorrectorRlx(A, b, C, y, S, M, dy1, muprimal, ncorr, delta, rhouser, bound);
    sl = y - l;
    su = u - y;
    
    ncorr = dsdpParam{27};
    
    % Logging
    nrmtk = abs(tau * kappa - muHSD);
    dObj = b' * y / tau;
    fprintf("%3d  %10.2e  %10.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %-8s\n",...
        i, pObj, dObj, nrmtk, muprimal, pinfeas, step, delta, newub);
    
    if pObj < dObj && isfirst
        pObj = 1e+10 + dObj;
        isfirst = false;
        if dsdpcheckDualRay(A, b, dymuprimal)
            reason = "DSDP_PRIMAL_INFEASIBLE_DUAL_UNBOUNDED";
            showdash(ndash);
            fprintf("Dual ray detected from iteration \n");
            showdash(ndash);
            break;
        end % End if    
    end % End if
    
    if (pObj - dObj) / (abs(pObj) + abs(dObj) + 1) < tol
        if pObj < dObj && pObj < 1e+10
            reason = "DSDP_NUMERICAL_ERROR";
        elseif pObj < dObj && pObj == 1e+10
            reason = "DSDP_LARGE_DOBJ";
        else
            reason = "DSDP_OPTIMAL";
        end % End if
        break;
    end % End if
    
    if dObj > pinfeasbound
        reason = "DSDP_LARGE_DOBJ";
        break;
    end % End if
    
end % End for

showdash(ndash);
if reason == "DSDP_OPTIMAL"
    fprintf("Reason: DSDP_OPTIMAL \n");
elseif reason == "DSDP_NUMERICAL_ERROR"
    fprintf("Reason: DSDP_NUMERICAL_ERROR \n");
elseif reason == "DSDP_SMALL_STEP"
    fprintf("Reason: DSDP_SMALL_STEP \n");
    fprintf("Gap: %10.5e \n", (pObj - dObj) / (abs(pObj) + abs(dObj) + 1));
elseif reason == "DSDP_LARGE_DOBJ"
    fprintf("Reason: DSDP_LARGE_DOBJ \n");
    fprintf("Suspecting dual unboundedness \n");
    if dsdpcheckDualRay(A, b, dymuprimal)
        reason = "DSDP_PRIMAL_INFEASIBLE_DUAL_UNBOUNDED";
        fprintf("Dual ray detected from iteration \n");
    end % End if
elseif reason == "DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE"
    fprintf("Reason: DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE \n");
    fprintf("Suspecting primal infeasibility \n");
    if dsdpcheckDualRay(A, b, dymuprimal)
        reason = "DSDP_PRIMAL_INFEASIBLE_DUAL_UNBOUNDED";
        fprintf("Dual ray detected from iteration \n");
    end % End if
end % End if
fprintf("Phase 2 finishes \n");
showdash(ndash);
end % End function