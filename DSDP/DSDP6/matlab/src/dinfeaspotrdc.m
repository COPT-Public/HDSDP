function [y, S] = dinfeaspotrdc(A, b, C, y, S, Ry, M, dy1, mu, ncorr)
% Potential reduction with presence of dual infeasibility

m = length(y);
dy2 = zeros(m, 1);
bz = zeros(m, 1);
for iter = 1:ncorr
    
    [L, ~, pmt] = lchol(S);
    pinv = dsdpInvPerm(pmt);
    
    for i = 1:m
        ASinv = cholSolve(L, pmt, pinv, A{i})';
        dy2(i) = trace(ASinv);
    end % End for
    
    sl = y + 1e+07;
    su = 1e+07 - y;
    
    % dy2 = dy2 - sl.^-1 + su.^-1;
    dy2 = M \ dy2;
    
    dycorr = 0 * dy1 / mu - dy2;
    dsl = dycorr;
    dsu = - dycorr;
    
    dScorr = - dsdpgetATy(A, dycorr);
    alpha = dsdpgetalpha(S, dScorr);
    
    if alpha < 0
        alpha = 1000;
    else
        alpha = min(0.95 * alpha, 1.0); 
    end
    
    step = - 1 / min(dsl ./ sl);
    if step > 0
        alpha = min(alpha, 0.95 * step);
    end % End if
    
    step = - 1 / min(dsu ./ su);
    if step > 0
        alpha = min(alpha, 0.95 * step);
    end % End if
    
    alpha = min(alpha, 1.0);
    
    oldpot = dsdpgetMeritValue(bz, y, mu, L);
    % oldpot = dsdpgetMeritValueRlx(bz, y, mu, L);
    btdy = 0; % b' * dycorr;
    
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
        % printf(pot);
        % pot = dsdpgetMeritValueRlx(bz, ynew, mu, L);
        
        anum = 2 * (pot - oldpot + btdy * alpha) / (alpha^2);
        bnum = btdy;
        
        if pot - oldpot > - abs(0.1 * alpha * btdy)
            if bnum / anum < alpha && bnum / anum > 0
                alpha = bnum / anum;
            else
                alpha = alpha * 0.5;
            end % End if
        else
            y = ynew;
            S = Snew;
            break;
        end % End if 
        
    end % End while
    
    if alpha <= 1e-04
        break;
    end % End if 
    
end % End for

end % End function