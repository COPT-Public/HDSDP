function [potred, pot] = getpotreduce(rho, A, ATA, x, potold, coneidx)

[f, ~] = fpot(A, ATA, x);
pot = rho * log(f) - sum(log(x(coneidx)));
potred = pot - potold;

end % End function