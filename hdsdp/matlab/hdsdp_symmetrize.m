function [S] = hdsdp_symmetrize(A)

dA = diag(A);
S = A + A' - diag(dA);

end % End function