function [x] = potRecur(A, ydim, maxiter, P, pc, pidx, ridx, Aorig)
% Implement potential reduction for LPs
% The first ydim columns in A refer to y

warning off;
[~, n] = size(A); m = ydim;
ncone = n - m; coneidx = ydim + 1 : n;
projidx = coneidx;
% projidx = n - 1 : n;
nproj = length(projidx);
rho = 1.1 * (ncone + sqrt(ncone));
e = ones(nproj, 1);

rng(24);
x_prev = zeros(n, 1);
x_prev(coneidx) = 1.0;

Abackup = A;
ATAbackup = A' * A;
AT = A';
ATA = ATAbackup;
PPT = P * P';
% One step potential reduction
[f, ~] = fpot(A, ATA, x_prev);
potold = rho * log(f) - sum(log(x_prev(coneidx)));
x_pres = x_prev;
recompute = true;

fvals = zeros(maxiter, 1);
potrds = zeros(maxiter, 1);

betamax = 1.0;
beta = betamax;

tic;
fprintf("%5s  %8s  %8s  %10s  %10s %10s %10s\n",...
        "iter", "fval", "pot", "alphag", "alpham", "beta", "potred");

x_cum = x_pres;

logstar = "";
allowcurv = 0;
usecurvature = 1;
curvinterval = 0;

alpha = [1, 1];
freq = 64;
Xbuff = zeros(n, freq);
rcount = 0;

potred = -100;

for i = 1:maxiter
    
    rcount = rcount + 1;
    Xbuff(:, rcount) = x_pres;
%     assert(min(x_pres(coneidx)) > 0);

    % Anderson acceleration
    if mod(rcount, freq) == 0
        Xbuff(1:m, 4) = 0;
        Xbuff(coneidx, 4) = 1;
        Xbuff(1:m, 1:3) = randn(m, 3);
        Xbuff(coneidx, 1:3) = rand(ncone, 3);
        Xbuff(coneidx, 1) = Xbuff(coneidx, 1) / sum(Xbuff(coneidx, 1)) * ncone;
        Xbuff(coneidx, 2) = Xbuff(coneidx, 2) / sum(Xbuff(coneidx, 2)) * ncone;
        Xbuff(coneidx, 3) = Xbuff(coneidx, 3) / sum(Xbuff(coneidx, 3)) * ncone;
        
        AX = A * Xbuff;
        Q = AX' * AX;
        model.Q = sparse(Q) + speye(64) * 1e-08;
        model.A = sparse([Xbuff(coneidx, :); 
                          ones(1, freq)]);
        model.rhs = zeros(ncone + 1, 1);
        model.rhs(end) = 1.0;
        model.sense = [repelem('>', ncone, 1); '='];
        model.lb = -inf(freq, 1);
        param.BarHomogeneous = 1;
        param.OutputFlag = 0;
        grbsol = gurobi(model, param);
        andalp = grbsol.x;
%         cvx_begin quiet
%         cvx_solver gurobi
%         cvx_precision(1e-10)
%         variable andalp(freq, 1)
%         minimize( quad_form(andalp, AX' * AX) );
%         subject to
%             sum(andalp) == 1.0;
%             Xbuff(coneidx, :) * andalp >= 1e-10;
%         cvx_end

        xand = Xbuff * andalp;
        potold = 1e+10;
        fv = norm(A * x_pres);
        fvand = norm(A * xand);  
        
        xand(coneidx) = max(xand(coneidx), 1e-12);
        x_pres = xand;
        x_cum = x_cum .* x_pres;
        
        A = A * diag([ones(ydim, 1); x_pres(coneidx)]);
        AT = A';
        ATA = A' * A;
        
        x_prev(coneidx) = x_prev(coneidx) ./ x_pres(coneidx);
        x_pres(coneidx) = 1.0;
        
        rcount = 0;
    end % End if
        
    if mod(i, 500) == 0 && 0
%         fprintf("Restart \n");
%         keyboard;
        % Projective transformation
%         x_pres = sum((1:64)' .* Xbuff)' / ((rcount * (rcount + 1)) / 2);
%         x_prev = x_pres;
%         rho = rho * 1.5;
        potold = 1e+10;
        x_cum = x_cum .* x_pres;
        A = A * diag([ones(ydim, 1); x_pres(coneidx)]);
        AT = A';
        ATA = A' * A;
        x_prev(coneidx) = x_prev(coneidx) ./ x_pres(coneidx);
        x_pres(coneidx) = 1.0;
        rcount = 0;
    end % End if
    
    % Start potential reduction
    if recompute
        
        [f, g] = fpot(A, ATA, x_pres);
        
        fvals(i) = f;
        
        % Prepare momentum
        mk = x_pres - x_prev;
        % mk = (randn(n, 1) * 1e-02) .* mk + mk;
        
        % Gradient projection
        gk = g;
        gk(coneidx) = gk(coneidx) - (f ./ x_pres(coneidx)) / rho;
        
        if usecurvature && allowcurv
            logstar = "*";
            method = "direct";
            [mk, ~] = findnegacurv(x_pres, m, coneidx, projidx, rho, g, f, ATA, AT, A, [], method);
            usecurvature = false;
        end % End if 
        
        gk(projidx) = gk(projidx) - e * sum(gk(projidx)) / nproj;
        mk(projidx) = mk(projidx) - e * sum(mk(projidx)) / nproj;
        
        % Prepare hessian
        nrmgk = norm(gk);
        nrmmk = norm(mk);
        
        if nrmmk > 0.0
            mk = mk / nrmmk;
        end % End if
        
        gk = gk / nrmgk;
        
        nrmgk = nrmgk * rho / f;
        
        Agk = A * gk; Amk = A * mk;
        
        gTgk = g' * gk; gTmk = g' * mk;
        xinvgk = gk ./ x_pres; xinvmk = mk ./ x_pres;
        xinvgk(1:m) = 0; xinvmk(1:m) = 0;
        
        gkXXgk = norm(xinvgk)^2;
        mkXXmk = norm(xinvmk)^2;
        gkXXmk = xinvgk' * xinvmk;
        
        gkHgk = -rho * (gTgk / f)^2 + gkXXgk + rho * norm(Agk)^2 / f;
        mkHmk = -rho * (gTmk / f)^2 + mkXXmk + rho * norm(Amk)^2 / f;
        mkHgk = -rho * (gTmk * gTgk) / f^2 + gkXXmk + rho * Amk' * Agk / f;
        
        gkTgk = gk' * gk * nrmgk; gkTmk = gk' * mk * nrmgk;
        
        H = [gkHgk, mkHgk;
            mkHgk, mkHmk];
        h = [gkTgk; gkTmk];
        M = [gkXXgk, gkXXmk;
            gkXXmk, mkXXmk];
        
    end % End if
    
    [alpha, mval] = subtrust(H, h, M, beta^2 / 2.25, 1e-10);
    d = alpha(1) * gk + alpha(2) * mk;
    
    xtmp = x_pres + d;
    [ftmp, ~] = fpot(A, ATA, xtmp);
    potnew = rho * log(ftmp) - sum(log(xtmp(coneidx)));
    potred = potnew - potold;
    
    potrds(i) = potred;
    
    ratio = potred / mval;
    
    if isnan(ratio)
        ratio = 0;
    end % End if
    
    if f < 1e-20 || beta < 1e-05
        break;
    end % End if
    
    if (ratio < 0.2 || potred > 0)
        beta = beta * 0.25;
        recompute = false;
        continue;
    elseif ratio > 0.75
        beta = min(beta * 2, betamax);
    else
        
    end % End if
    
    recompute = true;
    x_prev = x_pres;
    x_pres = x_pres + d;
    potold = potnew;
    
    if mod(i, 500)
        fprintf("%5d  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e %1.1s %3.3f %+3.3e\n",...
            i, f, potnew, alpha(1), alpha(2), beta, potred, logstar, toc, min(x_pres(coneidx)));
    end % End if
    
    logstar = "";
    
    curvinterval = curvinterval + 1;
    if potred > -0.1 && curvinterval > 100
        usecurvature = true;
        curvinterval = 0;
    end % End if
    
%     x_pres(1:m) = PPT \ P * (-x_pres(pidx) + pc * x_pres(end));
    
end % End for

fprintf("%5d  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %3.3f\n",...
    i, f, potnew, alpha(1), alpha(2), beta, potred, toc);

x_pres(coneidx) = x_pres(coneidx) .* x_cum(coneidx);
x = x_pres;
% semilogy(fvals, 'LineWidth', 3);
% hold on;

end % End function