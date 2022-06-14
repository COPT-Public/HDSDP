function [y, S, pObj, muprimal, xmaker, reason] = dsdpPrimalFeasPhase(A, b, C, y, S, muHSD, pObj, dObj, tau, kappa, ispdfeas, dsdpParam)
% Implementing primal feasibility certificate for DSDP
% y and S are input such that ATy + S - C * tau = 0
% As if using original DSDP to solve
% maximize tau * b' * y s.t. ATy + S - C * tau = 0
[n, ~] = size(S);
m = length(y);
rhouser      = dsdpParam{6};
maxiter      = dsdpParam{1};
maxpfeasiter = dsdpParam{17};
ncorr        = dsdpParam{27};
tol          = dsdpParam{2};
ndash        = dsdpParam{20};
muprimal     = min((pObj - dObj) / rhouser, muHSD);
tau          = 1;
normalizer   = dsdpParam{23};
pinfeasbound = dsdpParam{26};
dObj = b' * y;
isfirst      = true;

if ispdfeas
    maxpfeasiter = dsdpParam{1};
end % End if

xmaker = cell(3, 1);
smallstep = zeros(3, 1);

reason = "DSDP_MAXITER";

for i = 1:maxiter
    
    if i >= maxpfeasiter + 1 && reason == "DSDP_MAXITER"
        reason = "DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE";
        break;
    end % End if
    
    ismufeas = false;
    muHSD = min(muHSD, muprimal);
    
    [M, ~, asinv, ~, ~, ~, ~, ~, ~, ~] = dsdpgetSchur(A, S);
    
    diagM = max(diag(M)); 
    Mhat = M;
    
    if muprimal < 1e-05
        % Mhat = M + eye(m, m) * 1e-08;
    end % End if
    
    if isfirst
        Mhat = Mhat + eye(m, m) * min(1e-08, diagM * 1e-06);
    end % End if
        
    % dy1 = dsdpConjGrad(Mhat, b, diagM, zeros(m, 1)) / tau;
    % dy2 = dsdpConjGrad(Mhat, asinv, diagM, zeros(m, 1));
    dy1dy2 = Mhat \ [b, asinv];
    dy1 = dy1dy2(:, 1) / tau;
    dy2 = dy1dy2(:, 2);
    
    dymuprimal = dy1 / muprimal - dy2;
    backwardnewton = C - dsdpgetATy(A, y - dymuprimal);
    
    % Proximity
    delta = sqrt(dymuprimal' * Mhat * dymuprimal);
    
    % Check feasibility of backward Newton step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if dsdpIspsd(backwardnewton)
        ismufeas = true;
        reason = "DSDP_PRIMAL_DUAL_FEASIBLE";
        muub = muprimal * (dymuprimal' * asinv + n + 2 * m);
        pObj = dObj + muub;
        muub = muub / n;
        mulb = muub / rhouser;
        newub = "*";
        xmaker{1} = y;
        xmaker{2} = dymuprimal;
        xmaker{3} = muprimal;
    else
        muub = (pObj - dObj) / (n + 2 * m);
        mulb = muub / rhouser;
        newub = "";
    end % End if
    
    % Select new mu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newmu = dsdpselectMu(A, S, muprimal, dymuprimal, dy1, backwardnewton, ismufeas, (pObj - dObj) / (n + 2 * m));
    muprimal = min(newmu, muub);
    muprimal = max(muprimal, mulb);
    
    if delta < 0.1
        muprimal = muprimal * 0.1;
    end % End if
    
    % Potential reduction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dymuprimal = dy1 / muprimal - dy2;
    dS = - dsdpgetATy(A, dymuprimal);
    [y, S, step] = dsdptakedualStep(A, b, C, pObj, y, dymuprimal, S, dS, (pObj - dObj) / muprimal);
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
    
    [y, S, muprimal] = dsdpdualCorrector(A, b, C, y, S, M, dy1, muprimal, ncorr, delta, rhouser);
    
    ncorr = dsdpParam{27};
    
    % Logging
    nrmtk = abs(tau * kappa - muHSD);
    dObj = b' * y / tau;
    fprintf("%3d  %10.2e  %10.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %-8s\n",...
        i, pObj, dObj, nrmtk, muprimal, muHSD, step, delta, newub);
    
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
    
    nrm = norm(y, 'inf');
    if norm(y, 'inf') > normalizer
        nrm = sqrt(nrm);
        y = y / nrm;
        S = S / nrm;
        C = C / nrm;
        tau = tau / nrm;
        fprintf("Norm y is large \n");
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