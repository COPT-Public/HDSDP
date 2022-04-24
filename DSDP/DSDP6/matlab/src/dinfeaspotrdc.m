function [y, S, Ry] = dinfeaspotrdc(A, b, C, y, S, Ry, M, dy1, mu, ncorr, bound, dratemax)
% Potential reduction with presence of dual infeasibility
% Here we allow the corrector to reduce dual infeasibility by incorporating
% the corresponding step as a component of the linesearch
if nargin < 11
    dratemax = 0.0;
end % End if

m = length(y);
dy2 = zeros(m, 1);
bz = zeros(m, 1);
d4 = bz;
asinvrysinv = zeros(m, 1);

if Ry(1, 1) == 0.0
    dratemax = 0.0;
    % ncorr = 1;
end % End if 

for iter = 1:ncorr
    
    [L, ~, pmt] = lchol(S);
    pinv = dsdpInvPerm(pmt);
    
    for i = 1:m
        ASinv = cholSolve(L, pmt, pinv, A{i})';
        dy2(i) = trace(ASinv);
        
        if dratemax
            SinvASinv = cholSolve(L, pmt, pinv, ASinv);
            asinvrysinv(i) = trace(SinvASinv) * Ry(1, 1);
        end % End if 
        
    end % End for
    
    sl = y + bound; su = bound - y;
    dy2 = dy2 - sl.^-1 + su.^-1;
    dy2 = M \ dy2;
    
    % An extra term of dual infeasibility
    if dratemax
        d4  = M \ asinvrysinv;
    end % End if 
    
    dycorr = - dy2;
    dsl = dycorr;
    dsu = - dycorr;
    
    dScorr = - dsdpgetATy(A, dycorr);
    alpha = dsdpgetalpha(S, dScorr);
    
    if alpha < 0
        alpha = 1000;
    else
        alpha = min(0.95 * alpha, 100); 
    end
    
    step = - 1 / min(dsl ./ sl);
    if step > 0
        alpha = min(alpha, step);
    end % End if
    
    step = - 1 / min(dsu ./ su);
    if step > 0
        alpha = min(alpha, step);
    end % End if
    
    alpha = min(0.95 * alpha, 1.0);
    
    oldpot = dsdpgetMeritValue(bz, y, mu, L);
    % oldpot = dsdpgetMeritValueRlx(bz, y, mu, L);
    
    while alpha > 1e-04
        
        ynew = y + alpha * dycorr;
        Snew = C - dsdpgetATy(A, ynew) - Ry;
        
        try
            L = lchol(Snew);
        catch
            alpha = alpha * 0.5;
            continue;
        end % End try
        
        pot = dsdpgetMeritValue(bz, ynew, mu, L);
        
        if pot - oldpot >= 0
           alpha = alpha * 0.5;
        else
            yhat = ynew; 
            Shat = Snew; 
            assert(min(y) > -1e+07);
            break;
        end % End if 
    end % End while
    
    if alpha <= 1e-04
        break;
    end % End if 
    
    if dratemax == 0.0
        y = yhat;
        S = Shat;
    else
        % alpha improves centrality
        gammafeas = dsdpgetalpha(Shat, Ry - dsdpgetATy(A, alpha * d4)); 
        drate = min(gammafeas * dratemax, 1.0);
        
        % Assemble the two directions and take step
        dycorr = - dy2 + drate * d4;
        y = y + alpha * dycorr;
        Ry = Ry * (1 - alpha * drate);
        S = -dsdpgetATy(A, y) + C - Ry;
        
        if Ry(1, 1) == 0.0
            break;
        end % End if
        
        if drate * alpha < 0.1
            dratemax = 0.0;
            ncorr = min(ncorr, iter + 2);
        end % End if
        
        % assert(dsdpIspsd(S));
        assert(max(y) < 1e+07);
        assert(min(y) > -1e+07);
    end % End if
    
end % End for

end % End function