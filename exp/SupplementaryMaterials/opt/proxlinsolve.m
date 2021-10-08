function [x] = proxlinsolve(a, b, xk, gamma, tau, maxiter)
% This function solves the proximal point sub-problem using prox-linear
% method

x_before = xk;
x_after = xk;

for i = 1:maxiter
    
    x_before = x_after;
    [batchtemp, n] = size(a);
    
    aTx = a * x_after;
    % Apply QP to obtain next iterate
    % There are batch + n variables
    qpn = batchtemp + n;
    
    % Construct quadratic term
    model.Q = speye(qpn) * (tau + gamma) / 2;
    model.Q(1:batchtemp, 1:batchtemp) = 0;
    
    % Construct linear coefficient
    model.obj = zeros(qpn, 1);
    model.obj(1:batchtemp) = 1 / batchtemp;
    model.obj(batchtemp + 1:end) = gamma * (x_after - xk);
    
    % Construct RHS
    model.rhs = zeros(2 * batchtemp, 1);
    model.rhs(1:batchtemp) = b - aTx.^2;
    model.rhs(batchtemp + 1:end) = model.rhs(1:batchtemp);
    
    % Get lower bound
    model.lb = - inf(qpn, 1);
    
    % Construct linear constraint matrix
    model.A = [speye(batchtemp), 2 * (aTx.*a);
        -speye(batchtemp), 2 * (aTx.*a)];
    
    % Get sign information
    model.sense = [repmat('>', batchtemp, 1); repmat('<', batchtemp, 1)];
    grbparam.OutputFlag = 0;
    grbparam.LogtoConsole = 0;
    sol = gurobi(model, grbparam);    
    x_after = sol.x(batchtemp + 1:end) + x_before;
    
    if norm(x_before - x_after) <= 1e-04
       break;
    end % End if
    
    % fprintf("Residual: %3.e \n", norm(x_before - x_after));
end % End for

x = x_after;
end % End function