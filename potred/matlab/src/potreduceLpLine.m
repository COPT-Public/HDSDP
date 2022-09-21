function [x] = potreduceLpLine(A, ydim, maxiter, ncurv, linesearch, printlevel)
% Implement potential reduction for LPs
% The first ydim columns in A refer to y

warning off;
[~, n] = size(A); m = ydim;
ncone = n - m; coneidx = ydim + 1 : n;
rho = (ncone + sqrt(ncone));
e = ones(ncone, 1);
x_prev = zeros(n, 1);
x_prev(coneidx) = e / ncone;
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
    fprintf("%5s  %8s  %8s  %10s  %10s %10s %10s\n", "iter", "fval", "pot", "alphag", "alpham", "beta", "potred");
end % End if

usecurvature = true;
logstar = "*";
nstar = 0;

for i = 1:maxiter
    
    % Start potential reduction
    if recompute
        
        [f, g] = fpot(A, ATA, x_pres);
        
        % Prepare momentum
        mk = x_pres - x_prev;
        
        % Gradient projection
        gk = rho / f * g;
        gk(coneidx) = gk(coneidx) - e ./ x_pres(coneidx);
        
        if usecurvature && ncurv
            nstar = nstar + 1;
            logstar = "*";
            d = zeros(n, 1);
            d(coneidx) = x_pres(coneidx).^-2;
            % Consider negative curvature of Hessian
            Hess = (- (g * g') / f + ATA) + diag(d * f / rho);
            H11 = Hess(1:m, 1:m); H12 = Hess(1:m, m+1:end);
            H22 = Hess(m + 1:end, m+1:end);
            PH12 = H12 - (H12 * e) * e' / ncone;
            HeeT = (H22 * e) * e' / ncone;
            PH22P = H22 - HeeT - HeeT' + sum(sum(H22)) / ncone^2;
            Hproj = [H11,   PH12;
                     PH12', PH22P'];
            % Hess = rho * (- (g * g') / f + ATA) + diag(f * d);
            [mk, emin] = eigs(Hproj, 1, 'smallestreal');
            if (mk' * Hproj * mk) > 0
%             if isnan(emin)
                [V, evals] = eig(Hproj, 'vector');
                [~, id] = min(evals);
                mk = V(:, id);
            end % End if
            usecurvature = false;
            
        end % End if 
        
        gk(coneidx) = gk(coneidx) - e * sum(gk(coneidx)) / ncone;
        mk(coneidx) = mk(coneidx) - e * sum(mk(coneidx)) / ncone;
        
        % Prepare hessian
        nrmgk = norm(gk);
        nrmmk = norm(mk);
        
        if nrmmk > 0.0
            mk = mk / nrmmk;
        end % End if
        
        gk = gk / nrmgk;
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
    
    while true
        [alpha, mval] = subtrust(H, h, M, beta^2 / 4, 1e-05);
        d = alpha(1) * gk + alpha(2) * mk;
        xtmp = x_pres + d;
        if min(xtmp(coneidx)) > 0
            break;
        end % End if
        beta = beta * 0.8;
    end % End while
    
    [ftmp, ~] = fpot(A, ATA, xtmp);
    potnew = rho * log(ftmp) - sum(log(xtmp(coneidx)));
    potred = potnew - potold;
    
    ratio = potred / mval;
    
    if isnan(ratio)
        ratio = 0;
    end % End if
    
    if f < 1e-16 || beta < 1e-05
        break;
    end % End if
    
    if (ratio < 0.2 || potred > 0)
        beta = beta * 0.25;
        recompute = false;
        continue;
    elseif ratio > 0.75
        beta = min(beta * 2, 10);
    else
        
    end % End if
    
    recompute = true;
    x_prev = x_pres;
    x_pres = x_pres + d;
    potold = potnew;
    
%     assert(min(x_pres(coneidx)) > 0);
%     assert((abs(sum(x_pres(coneidx)) - 1) < 1e-06));
    
    if logstar == "*" && potred > -0.5
        ncurv = false;
    end % End if 
    
    if mod(i, 500) && printlevel
        fprintf("%5d  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e %1.1s %3.3f\n",...
            i, f, potnew, alpha(1), alpha(2), beta, potred, logstar, toc);
    end % End if
    logstar = "";
    
    if potred > -0.1 || (mod(i, 100) == 1 && potred > -10)
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