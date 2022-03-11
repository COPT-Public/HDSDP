function [y, S] = dinfeaspotrdc(A, b, C, y, S, Ry, M, dy1, mu)
% Potential reduction with presence of dual infeasibility

m = length(y);
ncorr = 4;

dy2 = zeros(m, 1);

for iter = 1:ncorr
    
    [L, ~, pmt] = lchol(S);
    pinv = dsdpInvPerm(pmt);
    
    for i = 1:m
        ASinv = cholSolve(L, pmt, pinv, A{i})';
        dy2(i) = trace(ASinv);
    end % End for
    
    dy2 = M \ dy2;
    
    dycorr = dy1 / mu - dy2;
    dScorr = - dsdpgetATy(A, dycorr);
    alpha = dsdpgetalpha(S, dScorr);
    alpha = min(0.95 * alpha, 1.0);
    
    % oldpot = dsdpgetMeritValue(b, y, mu, L);
    oldpot = dsdpgetMeritValue(0, 0, mu, L);
    btdy = b' * dycorr;
    
    while alpha > 1e-04
        
        ynew = y + alpha * dycorr;
        Snew = C - dsdpgetATy(A, ynew) - Ry;
        try
            L = lchol(Snew);
        catch
            alpha = alpha * 0.5;
            continue;
        end % End try
        % pot = dsdpgetMeritValue(b, ynew, mu, L);
        pot = dsdpgetMeritValue(0, 0, mu, L);
        
        anum = 2 * (pot - oldpot + btdy * alpha) / (alpha^2);
        bnum = btdy;
        
        if pot - oldpot > - 0.1 * alpha * btdy
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