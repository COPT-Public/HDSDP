function [delta, pObj] = dsdpgetProxMeasure(d2, d3, tau, mu, dObj, M, A, C, y, Ry, asinvrysinv, rysinv, asinv, S, b)
% Get measure of proximity in DSDP
if nargin >= 14
    verifypfeas = true;
else
    A = [];
    C = [];
    verifypfeas = false;
end % End if

m = length(y);
[n, ~] = size(C);

dy1 = d2; % M \ b
dy2 = d3; % M \ asinv

dydelta = dy1 * tau / mu - dy2;
delta = dydelta' * M * dydelta;

if verifypfeas
    bwnt = y - dydelta;
    if (max(bwnt) > 1e+07 && min(bwnt) < -1e+07)
        pObj = inf; 
        return;
    end % End if 
        
    if dsdpIspsd((- Ry + C * tau - dsdpgetATy(A, bwnt)))
        pObj = dObj + mu / tau * (rysinv + (asinvrysinv + asinv)' * dydelta + (n + 2 * m)); 
    else
        pObj = inf;
    end % End if 
end % End if 

end % End function