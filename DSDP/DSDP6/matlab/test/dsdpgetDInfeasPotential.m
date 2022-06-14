function [pot] = dsdpgetDInfeasPotential(rho, pObj, dObj, L, sl, su)
% Compute the infeasible dual potential function

pot = rho * log(pObj - dObj) - 2 * sum(log(diag(L))) - sum(log(sl)) - sum(log(su));

end % End function