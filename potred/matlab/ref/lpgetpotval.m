function [pot] = lpgetpotval(rho, f, yp, ym, x, s, kappa, tau)
% Evaluate the value of the potential function

pot = rho * log(f) - sum(log(yp)) - sum(log(ym)) - sum(log(x))...
        - sum(log(s)) - log(kappa) - log(tau);

end % End function