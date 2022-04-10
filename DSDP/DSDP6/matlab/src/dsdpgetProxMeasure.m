function [delta, pObj] = dsdpgetProxMeasure(d2, d3, tau, mu, dObj, M, A, C, y, Ry, asinvrysinv, rysinv, asinv, S, b, np, bound)
% Get measure of proximity in DSDP
verifypfeas = true;

m = length(y);
[n, ~] = size(C);

dy1 = d2; % M \ b
dy2 = d3; % M \ asinv

dydelta = dy1 * tau / mu - dy2;
delta = dydelta' * M * dydelta;

if verifypfeas
    bwnt = y - dydelta;
    if (max(bwnt) > bound && min(bwnt) < -bound)
        pObj = inf; 
        return;
    end % End if 
        
    if dsdpIspsd((- Ry + C * tau - dsdpgetATy(A, bwnt)))
        pObj = dObj + mu / tau * (rysinv + (asinvrysinv + asinv)' * dydelta + np); 
    else
        pObj = inf;
    end % End if 
end % End if 

end % End function