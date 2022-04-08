function [mu, pobj, pinfeas] = dsdpgetmualpha(asinv, b, csinv, ub, M, tol)
% Evaluate mu(alpha) by golden search
if nargin < 6
    tol = 1e-05;
end % End if 

gr = (sqrt(5) + 1) / 2;
mulb = 0;
muub = ub;
tol = min(tol, muub * tol) / gr;
diff = (muub - mulb) / gr;

while diff > tol
    
    c = muub - diff;
    d = mulb + diff;
    cobj = getgoldobj(c, csinv, b, asinv, M);
    dobj = getgoldobj(d, csinv, b, asinv, M);
    if cobj < dobj
        muub = d;
    else
        mulb = c;
    end % End if 
    
    diff = (muub - mulb) / gr;
    
    % fprintf("%e \n", diff);
end % End while

mu = (muub + mulb) / 2;
pobj = getgoldobj(mu, csinv, b, asinv, M);
pinfeas = norm(b - mu * asinv);

end % End function