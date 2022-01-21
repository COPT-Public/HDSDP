function [y, S, kappa, tau] = dsdpCorrector(A, b, C, y, S, kappa, tau, mu, M, corralpha, ncorr, stepstrategy)
% Implement the corrector step of DSDP. Usable with/without dual
% infeasibility

olddelta = inf;

for k = 1:ncorr

    [csinv, csinvcsinv, u, asinv] = dsdpgetCorrVec(A, C, S);
    d1_3 = M \ [b / mu + u, asinv];
    d1 = d1_3(:, 1);
    d3 = d1_3(:, 2);
    
    b1 = b - mu * u;
    taudenom =  b1' * d1 + mu * csinvcsinv + kappa / tau;
    
    if taudenom < 1e-10
        dtaucorr = 0.0;
    else
        dtaucorr = mu / tau - kappa + mu * csinv + b1' * d3;
        dtaucorr = dtaucorr / taudenom;
    end % End if
    
    dycorr = d1 * dtaucorr - d3;
    dScorr = C * dtaucorr - dsdpgetATy(A, dycorr);
    dkappacorr = - kappa + (mu / tau) - (kappa / tau) * dtaucorr;
    step = dsdpgetStepsize(S, dScorr, kappa, dkappacorr, tau, dtaucorr, stepstrategy, corralpha);
    
    y = y + step * dycorr;
    S = S + step * dScorr;
    kappa = kappa + step * dkappacorr;
    tau = tau + step * dtaucorr;
    
    delta = (dycorr' * M * dycorr) - 2 * u' * dycorr * dtaucorr + (csinvcsinv + (1 / tau^2)) * dtaucorr^2;
    if delta > olddelta || delta < 0.1
        break;
    end % End if 
    olddelta = delta;
end % End for 


end % End function