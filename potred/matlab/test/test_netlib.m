function [] = test_netlib(fname, maxiter, maxmn, minmn)

data = preprocess(fname);
A = data.A;
b = data.b;
c = data.c;

[m, n] = size(A);

% linesearch = false;
% neweigs = false;

if max(m, n) > maxmn || min(m, n) < minmn
   return;
end % End if

params.maxIter = maxiter;
params.maxRuizIter = 200;
params.maxTime = 600.0;
params.coefScal = 1;
[x, y, s] = potlp(A, b, c, params);

if abs(c' * x - b' * y) / (1 + abs(c' * x) + abs(b' * y)) > 1.1e-04
    params.coefScal = 0;
    [x, y, s] = potlp(A, b, c, params);
end % End if

% nrmb = norm(b, 1); nrmc = norm(c, 1);
% b = b / (nrmb + 1); c = c / (nrmc + 1);
% 
% HSDAA = [sparse(m, m), A, sparse(m, n), sparse(m, 1), -b;
%          -A',         sparse(n, n),  -speye(n), sparse(n, 1), c;
%          b',          -c',  sparse(1, n), -1, 0];
% 
% [D, E, HSDA] = ruizscale(HSDAA, 30);
% lpsol = potreduceLp(HSDA, m, maxiter, true, linesearch, neweigs, 0);
% sol = lpsol .* E;
% 
% kappa = sol(end - 1);
% tau = sol(end);
% y = sol(1:m);
% y = y / tau;
% s = sol(m + n + 1 : m + 2 * n) / tau;
% x = sol(m + 1: m + n) / tau;
% 
% fname = char(fname);
% 
% fprintf("%-12s %10.3e %10.3e  %10.3e  %10.3e \n", fname(11:end-8), c' * x, ...
%     norm(A * x - b) / (1 + norm(b, 1)), norm(A' * y + s - c) / (norm(c, 1) + 1), ...
%     (c' * x - b' * y) / (abs(c' * x) + abs(b' * y) + 1));

end % End function