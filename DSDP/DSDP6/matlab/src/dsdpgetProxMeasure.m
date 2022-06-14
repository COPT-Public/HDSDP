 function [delta, pObj] = dsdpgetProxMeasure(d2, d3, tau, mu, dObj, M, A, C, y, Ry, asinvrysinv, rysinv, asinv, S, b, np, bound)
% Get measure of proximity in DSDP
verifypfeas = true;

m = length(y);
[n, ~] = size(C);

dy1 = d2; % M \ b
dy2 = d3; % M \ asinv

dydelta = dy1 * tau / mu - dy2;
delta = dydelta' * (b * tau / mu - asinv);

if verifypfeas
    bwnt = y - dydelta;
    if (max(bwnt) > bound && min(bwnt) < -bound)
        pObj = inf; 
        return;
    end % End if 
    
    backwardnewton = - Ry + C * tau - dsdpgetATy(A, bwnt);
    
    if dsdpIspsd(backwardnewton)
        pObj = dObj + mu / tau * (rysinv + (asinvrysinv + asinv)' * dydelta + np); 
         if (true)
            backwardnewtonub = bound - (y - dydelta);
            backwardnewtonlb = bound + (y - dydelta);
            sl = y + bound;
            su = bound - y;
            Xtmp = mu * (S \ (S \ (backwardnewton))');
            xl = mu * (sl.^-2).*backwardnewtonlb;
            xu = mu * (su.^-2).*backwardnewtonub;
            ax = zeros(m, 1);
            for q = 1:m
                ax(q) = trace(A{q} * Xtmp);
            end % End for
            resi = ax - xl + xu - b;
            % assert(norm(resi) < 1e-03);
        end % End if
    else
        pObj = inf;
    end % End if 
end % End if 

end % End function