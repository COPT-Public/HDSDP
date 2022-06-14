function [M, u, asinv, L, p, CSinv, csinv, csinvcsinv, asinvrdsinv, csinvrdsinv, rdsinv] = dsdpgetSchur(A, S, C, Rd, initstrategy)
% Compute Schur components

m = length(A);
[n, ~] = size(S);
if nargin <= 3
    computeRd = false;
    asinvrdsinv = [];
else
    computeRd = true;
    asinvrdsinv = zeros(m, 1);
end % End if

computeCu = true;
if nargin == 2
    computeCu = false;
    CSinv = [];
    csinv = 0.0;
    csinvcsinv = 0.0;
end % End if

csinvrdsinv = 0.0;

% Prepare arrays
[M, asinv, u] = dsdpPrepareArray(m);

[L, ~, pmt] = lchol(S);
% [R, ~, pmt] = chol(S, 'vector');
pinv = dsdpInvPerm(pmt);
% L = R';

if computeCu
    CSinv = cholSolve(L, pmt, pinv, C)';
    csinv = trace(CSinv);
    csinvcsinv = trace(CSinv * CSinv);
end % End if

if computeRd
    rdsinv = trace(cholSolve(L, pmt, pinv, speye(n))) * Rd(1, 1);
    SinvCSinv = cholSolve(L, pmt, pinv, CSinv);
    if initstrategy == "IdS"
        csinvrdsinv = trace(SinvCSinv * Rd);
    else
        csinvrdsinv = trace(SinvCSinv) * Rd(1, 1);
    end % End if
end % End if

% Assemble M
for p = 1:m
    ASinv = cholSolve(L, pmt, pinv, A{p})';
    asinv(p) = trace(ASinv);
    SinvASinv = cholSolve(L, pmt, pinv, ASinv);
    
    if computeRd
        if initstrategy == "IdS"
            asinvrdsinv(p) = trace(SinvASinv * Rd); %#ok<*AGROW>
        else
            asinvrdsinv(p) = Rd(1, 1) * trace(SinvASinv);
        end % End if
    end % End if
    
    if computeCu
        u(p) = trace(SinvASinv * C);
    end % End if
    
    for q = 1:p
        M(p, q) = trace(SinvASinv * A{q});
        M(q, p) = M(p, q);
    end % End for
    
end % End for

end % End function