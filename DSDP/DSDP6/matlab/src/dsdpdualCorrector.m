function [y, S, mu] = dsdpdualCorrector(A, b, C, y, S, M, dy1, mu, ncorr, delta, rho)
% Implement the dual corrector for DSDP

if ncorr == 0
    return;
end % End if 

m = length(y);

[n, ~] = size(S);
dy2 = zeros(m, 1);
bTdy1 = b' * dy1;

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
        mu = min(mu, (bTdy1 / bTdy2));
    end % End if
    mu = mu * n / (n + sqrt(n));
    
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
    alphabest = alpha * 0.9;
    bestmerit = oldmerit;
    
    while alpha > 1e-06
        ynew = y + alpha * dycorr;
        Snew = C - dsdpgetATy(A, ynew);
        try
            L = lchol(Snew);
        catch 
            alpha = alpha * 0.5;
            alphabest = alpha;
            continue;
        end % End try
        
        newmerit = dsdpgetMeritValue(b, ynew, mu, L);
        
        if newmerit <= oldmerit - 0.1 * alpha * bTdycorr
            y = ynew;
            S = Snew;
            break;
        else
            if newmerit < bestmerit
                alphabest = alpha;
                bestmerit = newmerit;
            end % End if
        end % End if
        
        alpha = alpha * 0.5;
    end % End while
    
    if alpha < 1e-06
        if alphabest < 1e-06
            break;
        end % End if
        y = y + alphabest * dycorr;
        S = C - dsdpgetATy(A, y);
    end % End if
    
    if y' * M * y < 0.1
        break;
    end % End if 
    
end % End for

assert(dsdpIspsd(S));

end % End function
