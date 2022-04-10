function [y, S] = dinfeaspotrdc(A, b, C, y, S, Ry, M, dy1, mu, ncorr, bound)
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
    
    sl = y + bound;
    su = bound - y;
    
    dy2 = dy2 - sl.^-1 + su.^-1;
    dy2 = M \ dy2;
    
    if isnan(dy2(1))
        error("Nan here");
    end % End if 
    
    dycorr = - dy2;
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
        alpha = min(alpha, step);
    end % End if
    
    step = - 1 / min(dsu ./ su);
    if step > 0
        alpha = min(alpha, step);
    end % End if
    
    alpha = min(0.95 * alpha, 1.0);
    
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
        
        if pot - oldpot > - abs(0.1 * alpha * btdy)
           alpha = alpha * 0.5;
        else
            y = ynew;
            assert(min(y) > -1e+07);
            S = Snew;
            break;
        end % End if 
        
    end % End while
    
    if alpha <= 1e-04
        break;
    end % End if 
    
end % End for

end % End function