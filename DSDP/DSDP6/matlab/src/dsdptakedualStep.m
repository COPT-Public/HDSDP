function [y, S, step] = dsdptakedualStep(A, b, C, pObj, y, dy, S, dS, rho)
% Take dual Newton step by potential reduction
% Refer to DSDPYStepLineSearch for detailed implementation

step = dsdpgetalpha(S, dS, 0.95);
if step < 0
    alpha = 1;
else
    alpha = min(1, 0.95 * step);
end % End if 

assert(alpha > 0);

L = lchol(S);
oldpotential = dsdpgetPotentialValue(rho, pObj, b, y, L);
bestpotential = inf;

alphabest = 0.0;

while alpha > 1e-08
   ynew = y + alpha * dy;
   Snew = C - dsdpgetATy(A, ynew);
   % Snew = S + alpha * dS;
   L = dsdpIspsd(Snew);
   if ~ L
       alpha = alpha / 3;
       continue;
   end % End if 
   
   L = lchol(Snew);
   
   newpotential = dsdpgetPotentialValue(rho, pObj, b, ynew, L);
   if newpotential < oldpotential - 0.05
       y = ynew;
       S = Snew;
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
    % S = S + alphabest * dS;
end % End if 

end % End function