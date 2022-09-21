function [fval, grad] = fpot(A, ATA, x)

grad = ATA * x;
nrm = norm(A * x);
fval = 0.5 * nrm^2;

end % End function