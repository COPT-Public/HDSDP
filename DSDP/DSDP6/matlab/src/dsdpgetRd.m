function [Rd] = dsdpgetRd(A, y, S, C, tau)

Rd = - dsdpgetATy(A, y) - S + C * tau;

end % End function