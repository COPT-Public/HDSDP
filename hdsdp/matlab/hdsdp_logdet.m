function [logdet] = hdsdp_logdet(S)
% Compute log determinant of S using Cholesky

L = chol(S)';
logdet = full(2 * sum(log(diag(L))));

end % End function