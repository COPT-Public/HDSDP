function [y, S, mu] = dsdpdualCorrector(A, b, C, y, S, M, dy1, mu, ncorr, delta, rho)
% Implement the dual corrector for DSDP

if ncorr == 0
    return;
end % End if 

m = length(y);

[n, ~] = size(S);
dy2 = zeros(m, 1);
bTdy1 = b' * dy1;

muold = mu;

for k = 1:ncorr
    [L, ~, pmt] = lchol(S);
    pinv = dsdpInvPerm(pmt);
    
    for i = 1:m
        ASinv = cholSolve(L, pmt, pinv, A{i})';
        dy2(i) = trace(ASinv);
    end % End for
    
    dy2 = M \ dy2;
    bTdy2 = b' * dy2;
    
    if bTdy1 > 0 && bTdy2 > 0
        mu = min(mu, (bTdy1 / bTdy2) / mu * muold);
        mu = mu * n / (n + sqrt(n));
    end % End if
    
    oldmerit = dsdpgetMeritValue(b, y, mu, L);
    
    dycorr = dy1 / mu - dy2;
    bTdycorr = b' * dycorr;
    
    dScorr = - dsdpgetATy(A, dycorr);
    step = dsdpgetalpha(S, dScorr);
    
    if step < 0
        alpha = 1.0;
    else
        alpha = min(1, 0.95 * step);
    end % End if
    
    alpha = min(alpha, rho / delta);
    
    while alpha > 1e-04
        ynew = y + alpha * dycorr;
        Snew = C - dsdpgetATy(A, ynew);
        try
            L = lchol(Snew);
        catch 
            alpha = alpha * 0.5;
            continue;
        end % End try
        
        newmerit = dsdpgetMeritValue(b, ynew, mu, L);
        
        anum = 2 * (newmerit - oldmerit + bTdycorr * alpha) / (alpha^2);
        bnum = bTdycorr;
        
        if newmerit <= oldmerit - max(0.1 * alpha * bTdycorr, 0.01 * alpha)
            y = ynew;
            S = Snew;
            break;
        else
            if bnum / anum < alpha && bnum / anum > 0
                alpha = bnum / anum;
            else
                alpha = alpha * 0.5;
            end % End if
        end % End if
        
    end % End while
    
    if alpha <= 1e-04
        break;
    end % End if
    
end % End for

end % End function
