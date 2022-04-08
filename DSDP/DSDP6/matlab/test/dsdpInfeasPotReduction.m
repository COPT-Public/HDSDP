function [y, S, step] = dsdpInfeasPotReduction(A, b, C, y, dy, S, alphamax, rho, pObj, Rd, boundy, dpenalty, fixmu)
% Implement potential reduction for DSDP dual infeasibility phase
if nargin <= 12
    fixmu = 0;
end % End if 

m = length(y);
asinv = zeros(m, 1);

rd = - Rd(1, 1);
sl = y + boundy;
su = boundy - y;
L = chol(S);
oldpotential = rho * log(pObj - b' * y + rd * dpenalty);
oldpotential = oldpotential - 2 * sum(log(diag(L))) - sum(log(su)) - sum(log(sl));

alpha = alphamax;
binf = norm(b, 'inf');

while alpha > 1e-08
    
    ynew = y + alpha * dy;
    Snew = C - dsdpgetATy(A, ynew) - Rd * (1 - alpha);
    
    try 
        L = lchol(Snew);
    catch
        alpha = alpha * 0.3;
        continue;
    end % End try
    
    for i = 1:m
        asinv(i) = trace(L \ (L' \ A{i}));
    end % End for
    csinv = trace(L \ (L' \ C));
    
    if fixmu
        newpObj = sum(abs(asinv * fixmu - b)) * boundy + fixmu * csinv;
    else
        [mualpha, newpObj, pinfeas] = dsdpgetmualpha(asinv, b, csinv, binf / min(abs(asinv)), boundy, 1e-04);
    end % End if 
        
    sl = ynew + boundy;
    su = boundy - ynew;
    dObj = b' * ynew - rd * (1 - alpha) * dpenalty;
    newpotential = dsdpgetDInfeasPotential(rho, newpObj, dObj, L, sl, su);
    
    if newpotential < oldpotential - 0.05
        y = ynew;
        S = Snew;
        break;
    else
        alpha = alpha * 0.3;
    end % End if 

end % End while

step = alpha;

end % End function