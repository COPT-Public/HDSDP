function [D, E, Ascal] = l2scale(A)
% Implement the l2-scaling

Asqr = A.^2;
D = sqrt(sum(Asqr, 2));
E = sqrt(sum(Asqr, 1))';

Ascal = diag(D.^-1) * A * diag(E.^-1);

end % End function