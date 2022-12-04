function [x] = potreduceLp(A, ydim, maxiter, ncurv, linesearch, neweigs, printlevel, lpA, lpb, lpc)
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
x_prev = zeros(n, 1);
x_prev(coneidx) = 1;
ATA = A' * A;
AT = A';

% One step potential reduction
[f, ~] = fpot(A, ATA, x_prev);
potold = rho * log(f) - sum(log(x_prev(coneidx)));
x_pres = x_prev;
recompute = true;

beta = 1.0;

tic;
if printlevel
    fprintf("%5s  %8s  %8s  %10s  %10s %10s %10s\n",...
        "iter", "fval", "pot", "alphag", "alpham", "beta", "potred");
end % End if

usecurvature = true;
logstar = "*";
nstar = 0;
vk = [];
for i = 1:maxiter
    
    % Start potential reduction
    if recompute
        
        [f, g] = fpot(A, ATA, x_pres);
        
        % Prepare momentum
        mk = x_pres - x_prev;
        
        % Gradient projection
%         gk = rho / f * g;
%         gk(coneidx) = gk(coneidx) - e ./ x_pres(coneidx);

        gk = g;
        gk(coneidx) = gk(coneidx) - (f ./ x_pres(coneidx)) / rho;

        if (usecurvature && ncurv)
            nstar = nstar + 1;
            logstar = "*";
            % Consider negative curvature of Hessian
            % Hess = rho * (- (g * g') / f + ATA) + diag(f * d);
            method = "lanczos";
            [mk, vk] = findnegacurv(x_pres, m, coneidx, projidx, rho, g, f, ATA, AT, A, vk, method);
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
        
        if gkHgk < 0
%             fprintf("- \n");
        end % End if
        
        h = [gkTgk; gkTmk];
        
        M = [gkXXgk, gkXXmk;
             gkXXmk, mkXXmk];
        
    end % End if
    
    [alpha, mval] = subtrust(H, h, M, beta^2 / 4, 1e-05);
    d = alpha(1) * gk + alpha(2) * mk;
    
    xtmp = x_pres + d;
    [ftmp, ~] = fpot(A, ATA, xtmp);
    potnew = rho * log(ftmp) - sum(log(xtmp(coneidx)));
    potred = potnew - potold;
    
    % More aggressive step via line search
    if linesearch
        
        if logstar == "*"
            aggstep = - 0.9995 / min(d(coneidx) ./ x_pres(coneidx));
            decay = 0.9;
        else
            aggstep = - 0.8 / min(d(coneidx) ./ x_pres(coneidx));
            decay = 0.5;
        end % End if

        while aggstep > 1.0
            xtmp = x_pres + aggstep * d;
            [linereduce, pottmp] = getpotreduce(rho, A, ATA, xtmp, potold, coneidx);
            if linereduce < potred
                potred = linereduce;
                potnew = pottmp;
                d = aggstep * d;
                break;
            end % End if
            aggstep = aggstep * decay;
        end % End while
        
    end % End if
    
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
        beta = min(beta * 2, 0.99995);
    else
        
    end % End if
    
    recompute = true;
    x_prev = x_pres;
    x_pres = x_pres + d;
    potold = potnew;
    
    if logstar == "*" && potred > -0.1
        ncurv = false;
    end % End if 
    
    if mod(i, 500) && printlevel
        fprintf("%5d  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e %1.1s %3.3f %+3.3e\n",...
            i, f, potnew, alpha(1), alpha(2), beta, potred, logstar, toc, min(x_pres(coneidx)));
    end % End if
    logstar = "";
     
    if potred > -0.01 || (mod(i, 100) == 1 && potred > -10)
        usecurvature = true;
    end % End if
    
end % End for

if printlevel
    fprintf("%5d  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %3.3f\n",...
        i, f, potnew, alpha(1), alpha(2), beta, potred, toc);
    fprintf("Negative curvature is computed %d times \n", nstar);
end % End if
x = x_pres;

end % End function