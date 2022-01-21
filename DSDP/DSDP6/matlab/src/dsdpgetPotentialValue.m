function [potential] = dsdpgetPotentialValue(rho, pObj, b, y, L)
% Retrieve the value of the potential function 

potential = rho * log(pObj - b' * y) - 2 * sum(log(diag(L)));

end % End function