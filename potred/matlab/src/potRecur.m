function [x] = potRecur(A, ydim, maxiter)
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
x_prev(coneidx) = 1.0;

Abackup = A;
ATAbackup = A' * A;
AT = A';
ATA = ATAbackup;
% One step potential reduction
[f, ~] = fpot(A, ATA, x_prev);
potold = rho * log(f) - sum(log(x_prev(coneidx)));
x_pres = x_prev;
recompute = true;

beta = 1.0;

tic;
fprintf("%5s  %8s  %8s  %10s  %10s %10s %10s\n",...
        "iter", "fval", "pot", "alphag", "alpham", "beta", "potred");

x_cum = x_pres;

logstar = "";
allowcurv = false;
usecurvature = false;
curvinterval = 0;

alpha = [1, 1];
Xbuff = zeros(64, n);
rcount = 0;

for i = 1:maxiter
    
    rcount = rcount + 1;
%     Xbuff(rcount, :) = x_pres;
        
    if mod(i, 500) == 0
%         fprintf("Restart \n");
%         keyboard;
        % Projective transformation
%         x_pres = sum((1:64)' .* Xbuff)' / ((rcount * (rcount + 1)) / 2);
%         x_prev = x_pres;
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
        
        % Prepare momentum
        mk = x_pres - x_prev;
        
        % Gradient projection
        %         gk = rho / f * g;
        %         gk(coneidx) = gk(coneidx) - e ./ x_pres(coneidx);
        
        gk = g;
        gk(coneidx) = gk(coneidx) - (f ./ x_pres(coneidx)) / rho;
        
        if usecurvature && allowcurv
            logstar = "*";
            method = "lanczos";
            [mk, ~] = findnegacurv(x_pres, m, coneidx, projidx, rho, g, f, ATA, AT, A, [], method);
            usecurvature = false;
        end % End if 
         
%         vmin = min(x_pres(coneidx));
%         
%         if (vmin < 0.001 * f^0.5)
%             [vmin, idmin] = mink(x_pres(coneidx), 5);
%             [vmax, idmax] = maxk(x_pres(coneidx), 5);
%             mk = zeros(n, 1);
%             mk(idmin + m) = -1;
%             mk(idmax + m) = 1;
%         end % End if 
        
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
    
    [alpha, mval] = subtrust(H, h, M, beta^2 / 4, 1e-05);
    d = alpha(1) * gk + alpha(2) * mk;
    
    xtmp = x_pres + d;
    [ftmp, ~] = fpot(A, ATA, xtmp);
    potnew = rho * log(ftmp) - sum(log(xtmp(coneidx)));
    potred = potnew - potold;
    
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
    
    if mod(i, 500)
        fprintf("%5d  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e %1.1s %3.3f %+3.3e\n",...
            i, f, potnew, alpha(1), alpha(2), beta, potred, logstar, toc, min(x_pres(coneidx)));
    end % End if
    
    logstar = "";
    
    curvinterval = curvinterval + 1;
    if potred > -0.1 && curvinterval > 50
        usecurvature = true;
        curvinterval = 0;
    end % End if
    
end % End for

fprintf("%5d  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %3.3f\n",...
    i, f, potnew, alpha(1), alpha(2), beta, potred, toc);

x_pres(coneidx) = x_pres(coneidx) .* x_cum(coneidx);
x = x_pres;

end % End function