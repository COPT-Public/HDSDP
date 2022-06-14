%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Linear Programming             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

rng(24);

m = 10;
n = 500;
ntest = 100;

alpha = 0.9;
sigma = 0.9;
tol = 1e-07;

isgood = zeros(ntest, 1);

for i = 1:ntest

    A = [rand(m, n - m), eye(m)];
    b = rand(m, 1);
    c = rand(n, 1);

    [xbench, fval, ~, ~, lambda] = linprog(c, [], [], A, b, zeros(n, 1), []);


    [x, y, s, iter] = dLPKappaTauPds(A, b, c, alpha, sigma, tol, 100, false);
    % [x, y, s, iter] = dLPTauPds(A, b, c, alpha, sigma, tol, 100, false);
    % [x,y,s,tau,kappa,iter] = DualHOlptaukappaNoVerify(A,b,c);

    if abs(b' * y - fval) / fval < 1e-05
    % if norm(x - xbench) < 1e-06 
        isgood(i) = 1;
        fprintf("Case %d passed \n", i)
    else
        fprintf("Base case: %3.e \n", abs(b' * y - fval) / fval)
        % error("Bad case");
    end % End if

end % End for

fprintf("In all %d good\n", sum(isgood));