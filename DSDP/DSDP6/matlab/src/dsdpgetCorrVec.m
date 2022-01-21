function [csinv, csinvcsinv, u, asinv] = dsdpgetCorrVec(A, C, S)
% Compute u, d1, d3 and csinv for corrector step
m = length(A);

u = zeros(m, 1);
asinv = zeros(m, 1);
[R, ~, pmt] = chol(S, 'vector');
pinv = dsdpInvPerm(pmt);
L = R';

CSinv = cholSolve(L, pmt, pinv, C)';
csinv = trace(CSinv);
csinvcsinv = trace(CSinv * CSinv);
SinvCSinv = cholSolve(L, pmt, pinv, CSinv);

for i = 1:m    
    ASinv = cholSolve(L, pmt, pinv, A{i});
    u(i) = trace(A{i} * SinvCSinv);
    asinv(i) = trace(ASinv);
end % End for
    
end % End function