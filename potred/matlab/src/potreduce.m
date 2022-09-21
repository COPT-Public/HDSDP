function [x] = potreduce(A, maxiter, usescale)
% Implement potential reduction

[~, n] = size(A);

warning off;

rho = n + sqrt(n);
e = ones(n, 1);
x_prev = e / n;
ATA = A' * A;

% One step potential reduction
x = x_prev;
[f, g] = fpot(A, ATA, x);
potnew = rho * log(f) - sum(log(x));
gk = rho / f * g - e./x;
beta = 0.1; % 1 / (2 + rho * gamma / f);
lx = sum( (x.^2).* gk * f ) / ( norm(x)^2 * rho);
px = x.*(rho / f * (g - e * lx)) - e;
d = - beta / norm(px) * (x .* px);
x_pres = x_prev + d;
potold = potnew;
recompute = true;

beta = 1.0;
fprintf("%5s  %8s  %8s  %10s  %10s %10s %10s\n", "iter", "fval", "pot", "alphag", "alpham", "deltag", "potred");

for i = 1:maxiter
    % Start potential reduction
        
    if recompute
        [f, g] = fpot(A, ATA, x_pres);
        
        % Prepare momentum
        mk = x_pres - x_prev;
        
        % Prepare hessian
        if usescale
            lx = sum( (x_pres.^2).* gk * f ) / ( norm(x_pres)^2 * rho);
            px = x_pres.*(rho / f * (g - e * lx)) - e;
            gk = 1.0 / norm(px) * (x_pres .* px);
        else
            % Gradient projection
            gk = rho / f * g - e./x_pres;
            gk = gk - e * sum(gk) / n;
        end % End if
        
        nrmgk = norm(gk);
        nrmmk = norm(mk);
        
        mk = mk / nrmmk;
        gk = gk / nrmgk;

        Agk = A * gk; Amk = A * mk;
        
        gTgk = g' * gk; gTmk = g' * mk;
        xinvgk = gk ./ x_pres; xinvmk = mk ./ x_pres;
        
        gkXXgk = norm(xinvgk)^2;
        mkXXmk = norm(xinvmk)^2;
        gkXXmk = xinvgk' * xinvmk;
        
        gkHgk = -rho * (gTgk / f)^2 + gkXXgk + norm(Agk)^2 / f;
        mkHmk = -rho * (gTmk / f)^2 + mkXXmk + norm(Amk)^2 / f;
        mkHgk = -rho * (gTmk * gTgk) / f^2 + gkXXmk + Amk' * Agk / f;
        
        gkTgk = gk' * gk * nrmgk; gkTmk = gk' * mk * nrmgk;
        
        H = [gkHgk, mkHgk;
             mkHgk, mkHmk];
        
        h = [gkTgk; gkTmk];
        
        M = [gkXXgk, gkXXmk;
             gkXXmk, mkXXmk];
        
    end % End if
    
    [alpha, mval] = subtrust(H, h, M, beta^2 / 4, 1e-10);
    d = alpha(1) * gk + alpha(2) * mk;
    
    xtmp = x_pres + d;
    [ftmp, ~] = fpot(A, ATA, xtmp);
    potnew = rho * log(ftmp) - sum(log(xtmp));
    potred = potnew - potold;
    
    ratio = potred / mval;
    
    if isnan(ratio)
        ratio = 0;
    end % End if
    
    if f < 1e-20 || beta < 1e-04
        break;
    end % End if
    
    if ratio < 0.1
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
    
    assert(min(x_pres) > 0);
    
    if mod(i, 500)
        fprintf("%5d  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e\n",...
            i, potnew, f, alpha(1), alpha(2), nrmgk * alpha(1), potred);
    end % End if 
    
    
end % End for

x = x_pres;

end % End function