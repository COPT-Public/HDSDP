function [y, S, sl, su, step] = dsdptakedualStepRlx(A, b, C, l, u, pObj, y, dy, S, dS, sl, dsl, su, dsu, rho)
% Take dual Newton step by potential reduction
% Refer to DSDPYStepLineSearch for detailed implementation

step = dsdpgetalpha(S, dS);

if step < 0
    alpha = 1;
else
    alpha = min(1, 0.95 * step);
end % End if 

step = - 1 / min(dsl ./ sl);
if step > 0
    alpha = min(1, 0.95 * step);
end % End if 

step = - 1 / min(dsu ./ su);
if step > 0
    alpha = min(1, 0.95 * step);
end % End if 

assert(alpha > 0);

L = lchol(S);
oldpotential = dsdpgetPotentialValue(rho, pObj, b, y, L);
oldpotential = oldpotential - sum(log(sl)) - sum(log(su));
bestpotential = inf;

alphabest = 0.0;

while alpha > 1e-08
   ynew = y + alpha * dy;
   Snew = C - dsdpgetATy(A, ynew);
   sunew = u - ynew;
   slnew = -l + ynew;
   % Snew = S + alpha * dS;
   L = dsdpIspsd(Snew);
   if ~ L || min(sunew) < 0.0 || min(slnew) < 0.0
       alpha = alpha / 3;
       continue;
   end % End if 
   
   L = lchol(Snew);
   
   newpotential = dsdpgetPotentialValue(rho, pObj, b, ynew, L);
   newpotential = newpotential - sum(log(slnew)) - sum(log(sunew));
   
   if newpotential < oldpotential - 0.05
       y = ynew;
       S = Snew;
       sl = slnew;
       su = sunew;
       break;
   else
       if newpotential < bestpotential
           alphabest = alpha;
           bestpotential = newpotential;
       end % End if 
   end % End if 
   
   alpha = alpha * 0.3;
   
end % End while

step = alpha;

if step <= 1e-08
    step = alphabest;
    y = y + alphabest * dy;
    S = C - dsdpgetATy(A, y);
    sl = - l + y;
    su = u - y;
    % S = S + alphabest * dS;
end % End if 

end % End function