function [delta, pObj] = dsdpgetProxMeasure(d12, d2, d3, u, csinvcsinv, csinv, tau, mu, dObj, M, A, C, y, Ry)
% Get measure of proximity in DSDP
if nargin >= 14
    verifypfeas = true;
else
    A = [];
    C = [];
    verifypfeas = false;
end % End if

taudenom = - u' * d12 + csinvcsinv + (1 / tau^2);

if abs(taudenom) < 1e-15
    dtaudelta = 0.0;
else
    dtaudelta = - dObj + mu / tau + mu * csinv + tau * u' * d2 - mu * u' * d3;
    dtaudelta = (1 / mu) *  dtaudelta / taudenom;
end % End if 

dydelta = d12 * dtaudelta + d2 * tau / mu - d3;
delta = (dydelta' * M * dydelta) - 2 * u' * dydelta * dtaudelta + (csinvcsinv + (1 / tau^2)) * dtaudelta^2;
delta = sqrt(delta);

if verifypfeas
    if dsdpIspsd(- Ry + C * (tau - dtaudelta) - dsdpgetATy(A, y - dydelta))
        kappadelta = mu * tau^-2 * (tau - dtaudelta);
        pObj = dObj - kappadelta;
    else
        pObj = inf;
    end % End if 
end % End if 

end % End function