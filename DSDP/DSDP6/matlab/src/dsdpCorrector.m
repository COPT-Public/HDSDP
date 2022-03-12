function [y, S, kappa, tau] = dsdpCorrector(A, b, C, y, S, Ry, kappa, tau, mu, M, d2)
% Implement the corrector step of DSDP. Usable with/without dual
% infeasibility

ncorr = 4;
d11 = d2 / mu;

for k = 1:ncorr

    [L, ~, pmt] = lchol(S);
    
    [csinv, csinvcsinv, u, asinv] = dsdpgetCorrVec(A, C, L, pmt);
    
    d12_3 = M \ [u, asinv];
    d12 = d12_3(:, 1);
    d3 = d12_3(:, 2);
    
    taudenom = (b - mu * u)' * (d11 + d12) + mu * csinvcsinv + kappa / tau;
    
    if abs(taudenom) < 1e-20
        dtaucorr = 0.0;
    else
        dtaucorr = -b' * y + mu / tau + mu * csinv;
        dtaucorr = dtaucorr - (b' * d2) * (tau / mu) + b' * d3 + u' * d2 * tau;
        dtaucorr = dtaucorr - mu * u' * d3;
        dtaucorr = dtaucorr / taudenom;
    end % End if 
    
    dycorr = (d11 + d12) * dtaucorr + d2 * tau / mu - d3;
    dScorr = C * dtaucorr - dsdpgetATy(A, dycorr);
    dkappacorr = - kappa + (mu / tau) - (kappa / tau) * dtaucorr;
    
    bTdycorr = b' * dycorr * tau;
        
    step = dsdpgetStepsize(S, dScorr, kappa, dkappacorr, tau, dtaucorr, "xxx", 0.97);
    
    % Get potential - log det (S) - log tau - log kappa
    oldpot = dsdpgetMeritValue(0, 0, 1.0, L) - log(tau); %  - log(kappa);
    oldpot = oldpot * mu;
    
    while step > 1e-03
        
        ynew = y + step * dycorr;
        taunew = tau + step * dtaucorr;
        kappanew = kappa + step * dkappacorr;
        Snew = C * taunew - dsdpgetATy(A, ynew) - Ry;
    
        try 
            L = lchol(Snew);
        catch
            step = step * 0.5;
            continue;
        end % End try
        
        pot = dsdpgetMeritValue(0, 0, 1.0, L) - log(taunew); %  - log(kappanew);
        pot = pot * mu;
        
        anum = 2 * (pot - oldpot + bTdycorr * step);
        bnum = bTdycorr;
        
        if pot - oldpot > - 0.1 * step * bTdycorr
            if bnum / anum < step && bnum / anum > 0
                step = bnum / anum;
            else
                step = step * 0.5;
            end % End if
        else
            y = ynew;
            S = Snew;
            tau = taunew;
            kappa = kappanew;
            break;
        end % End if
    
    
        if step <= 1e-03
            break;
        end % End if 
        
    end % End while
        
end % End for 

end % End function