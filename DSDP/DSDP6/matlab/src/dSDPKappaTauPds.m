function [X, S, y, kappa, tau, reason] = dSDPKappaTauPds(A, b, C, dsdpParam)
% DSDP technique applied to SDP with HSD

maxiter      = dsdpParam{1};
tol          = dsdpParam{2};
sigma        = dsdpParam{3};
alpha        = dsdpParam{4};
rho          = dsdpParam{6};
initstrategy = dsdpParam{7};
reuse        = dsdpParam{8};
corr         = dsdpParam{9};
maxpfeasiter = dsdpParam{17};
initbeta     = dsdpParam{10};
ndash        = dsdpParam{20};
presolve     = dsdpParam{21};

iter = 999;
dinfeas = ["DSDP_DUAL_INFEASIBLE", "DSDP_PRIMAL_FEASIBLE_DUAL_INFEASIBLE"];

if ~ issparse(C)
    C = sparse(C);
end % End if

% Extract size
m = length(b);
[n, ~] = size(C);

% Prepare iteration arrays
mu      = dsdpParam{11};
tau     = dsdpParam{12};
kappa   = dsdpParam{13};

Rdflag = true;
reason = "DSDP_PRIMAL_DUAL_UNKNOWN";

Aold = A;
bold = b;
Cold = C;

if presolve
    fprintf("Start presolving \n");
    [A, b, C, pscaler, dscaler] = dsdpPresolver(A, b, C);
    fprintf("Presolving done \n");
end % End if

% Prepare initial solutions
[Rd, S, y] = dsdpInitialize(A, C, tau, initstrategy, initbeta);

nora = 0.0;
% if initstrategy == "IdS"
%     nrmRd = norm(Rd, 'fro');
% else
%     nrmRd = sqrt(n) * Rd(1, 1);
% end % End if

showdash(ndash);
fprintf("Phase 1 start. Eliminating dual infeasibility \n");
showdash(ndash);
fprintf("%4s  %10s  %10s  %8s  %8s  %8s  %8s  %8s  %8s \n", "Iter", "pObj", "dObj", "dInf", "tkInf", "muPrimal", "muHSD", "step", "pNorm");
showdash(ndash);

% Phase 1
[S, y, kappa, tau, muHSD, pObj, reason, iterphase1] =  dsdpDualInfeasPhase(A, b, C, y, S, Rd, mu, kappa, tau, dsdpParam);

fprintf("Reason: %s \n", reason);
showdash(ndash);

fprintf("Phase 1 finished. \n");
showdash(ndash);

if ismember(reason, dinfeas)
    S = S / kappa;
    y = y / kappa;
    X = sparse(n, n);
    iter = iterphase1;
    fprintf("DSDP Ends. Status: Dual Infeasibile \n");
    return
end % End if

% Certificate found for both primal and dual
if reason == "DSDP_PRIMAL_DUAL_FEASIBLE"
    fprintf("Phase 1 certificates primal and dual feasibility \n");
    S = S / tau;
    y = y / tau;
    dObj = b' * y;
    fprintf("Phase 2 start. Restarting dual scaling with pObj: %8.2e dObj: %8.2e \n", pObj, dObj);
    ispdfeas = true;
elseif reason == "DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE"
    fprintf("Phase 1 only certificates dual feasibility \n");
    Sheur = S / tau;
    yheur = y / tau;
    dObj = b' * yheur;
    pObj  = max(dsdpParam{5}, dObj + dsdpParam{5} / 10);
    fprintf("Phase 2 start. Trying to certificate primal feasibility in %d iterations \n", maxpfeasiter);
    ispdfeas = false;
end % End if

showdash(ndash);
fprintf("%4s  %10s  %10s  %8s  %8s  %8s  %8s  %8s \n", "Iter", "pObj", "dObj", "tkInf", "muPrimal", "muHSD", "step", "pNorm");
showdash(ndash);


% Phase 2
if ispdfeas
    [y, S, pObj, muprimal, xmaker, reason] = dsdpPrimalFeasPhase(A, b, C, y, S, muHSD, pObj, dObj, tau, kappa, ispdfeas, dsdpParam);
else
    [y, S, pObj, muprimal, xmaker, reason] = dsdpPrimalFeasPhase(A, b, C, yheur, Sheur, muHSD, pObj, dObj, tau, kappa, ispdfeas, dsdpParam);
end % End if

X = sparse(n, n);
if ismember(reason, ["DSDP_OPTIMAL", "DSDP_SMALL_STEP","DSDP_MAXITER"])
    
    if reason == "DSDP_MAXITER"
        fprintf("Maximum iteration reached \n");
    elseif reason == "DSDP_OPTIMAL"
        fprintf("DSDP Converges \n");
    end % End if 
    
    if ~ isempty(xmaker{1})
        ymaker  = xmaker{1};
        dymaker = xmaker{2};
        mumaker = xmaker{3};
        bnmaker = dsdpgetATy(A, dymaker);
        Smaker = C - dsdpgetATy(A, ymaker);
        R = chol(Smaker);
        D = R \ eye(n);
        % X = Smaker \ eye(n) + Smaker \ (Smaker \ bnmaker)';
        % X = X * mumaker;
        X = R \ (R \ (speye(n) + R' \ (R' \ bnmaker)'))';
        % X = D * (speye(n) + D' * bnmaker * D) * D';
        X = X * mumaker;
    end % End if 
    
end % End if

if reason == "DSDP_LARGE_DOBJ" || reason == "DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE"
    fprintf("Phase 2 suspects primal infeasibility \n");
    fprintf("Looking for dual improving ray \n");
    showdash(ndash);
    [y, S, reason] = dsdpprovepInfeas(A, b, C, y, S, kappa, 0, reason, dsdpParam);
    fprintf("Reason: %s \n", reason);
    if reason == "DSDP_PRIMAL_INFEASIBLE_DUAL_UNBOUNDED"
        fprintf("DSDP certificates primal infeasibility \n");
        showdash(ndash);
        fprintf("DSDP Ends \n");
        X = X / dscaler;
        return;
    else
        fprintf("Going back to Phase 2 \n");
    end % End if
end % End if

y = y ./ pscaler * dscaler;
S = S * dscaler;
    
pinfeas = zeros(m, 1);
Rd = dsdpgetATy(Aold, y) + S - Cold;

for i = 1:m
    pinfeas(i) = trace(Aold{i} * X) - bold(i);
end % End for

bn1 = norm(bold, 1);
cn1 = sum(sum(abs(Cold)));
err1 = norm(pinfeas) / (1 + bn1);
err2 = max(0, - eigs(X, 1, 'smallestreal')) / (1 + bn1);
err3 = norm(Rd, 'fro') / (1 + cn1);
err4 = max(0, - eigs(S, 1, 'smallestreal')) / (1 + cn1);
err5 = (trace(Cold * X) - bold' * y) / (1 + abs(trace(Cold * X)) + abs(bold' * y));
err6 = trace(X * S) / (1 + abs(trace(Cold * X)) + abs(bold' * y));

maxerr = max(abs([err1, err2, err3, err4, err5, err6]));
if maxerr > 1e-04 && ~ ismember(reason, [dinfeas, "DSDP_PRIMAL_INFEASIBLE_DUAL_UNBOUNDED"])
    showdash(ndash);
    fprintf("Inaccurate solution. DIMACS error exceeds %10.5e \n", maxerr);
    if max(abs(err5), abs(err6)) > 1e-02
        fprintf("Re-validating the solution \n");
        showdash(ndash);
        ynew = y .* pscaler / dscaler;
        Snew = C - dsdpgetATy(A, ynew);
        [~, ~, reason] = dsdpprovepInfeas(A, b, C, ynew, Snew, 100, 0, reason, dsdpParam);
        showdash(ndash);
        if reason == "DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE"
            fprintf("DSDP Failed \n");
        end % End if 
    end % End if 
end % End if 

fprintf("DSDP Ends \n");
fprintf("DIMACS Error: %e %e %e %e %e %e \n", err1, err2, err3, err4, err5, err6);
end % End function