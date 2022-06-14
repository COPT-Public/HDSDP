function [logdet] = dsdpgetlogdet(S)

L = lchol(S);
logdet = 2 * sum(log(diag(L)));

end % End function