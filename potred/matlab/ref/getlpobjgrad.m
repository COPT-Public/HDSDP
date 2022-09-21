function [f, g, pres, dres, compl] = getlpobjgrad(A, b, c, x, y, s, kappa, tau)
% Compute the function value, gradient and residuals

pres = A * x - b * tau;
dres = -A' * y - s + c * tau;
compl = b' * y - c' * x - kappa;

pnrm = norm(pres); dnrm = norm(dres); complnrm = abs(compl);
f = 0.5 * (pnrm^2 + dnrm^2 + complnrm^2);

g1 = -A * dres + b * compl;
g2 = A' * pres - c * compl;
g3 = -dres;
g4 = -compl;
g5 = -b' * pres + c' * dres;

g = [g1; g2; g3; g4; g5];      

end % End function