function [z] = proxlinblindsolve(u, v, b, zk, gamma, tau, maxiter)
% This function solves the proximal point sub-problem using prox-linear
% method

z_before = zk;
z_after = zk;
[batchtemp, n] = size(u);
xidx = 1:n;
yidx = (n + 1) : (2 * n);

for i = 1:maxiter
    
    z_before = z_after;
    uTx = u * z_after(xidx);
    vTy = v * z_after(yidx);
    
    gamma = gamma / alpha_0;
    % Get next iterate using Gurobi QP optimizer
    % Apply QP to obtain next iterate
    % There are batch + n variables
    qpn = batchtemp + 2 * n;
    
    % Construct quadratic term
    model.Q = speye(qpn) * (tau + gamma) / 2;
    model.Q(1:batchtemp, 1:batchtemp) = 0;
    
    % Construct linear coefficient
    model.obj = zeros(qpn, 1);
    model.obj(1:batchtemp) = 1 / batchtemp;
    
    % Construct RHS
    model.rhs = zeros(2 * batchtemp, 1);
    model.rhs(1:batchtemp) = b(batchidx) - (uTx .* vTy) - ...
        vTy .* (u * (w(xidx) - z_before(xidx))) - ...
        uTx .* (v * (w(yidx) - z_before(yidx)));
    model.rhs(batchtemp + 1:end) = - model.rhs(1:batchtemp);
    
    % Get lower bound
    model.lb = - inf(qpn, 1);
    
    % Construct linear constraint matrix
    model.A = [-speye(batchtemp), vTy.*u,  uTx.*v;
        -speye(batchtemp), - vTy.*u,  - uTx.*v];
    
    % Get sign information
    model.sense = repmat('<', 2 * batchtemp, 1);
    grbparam.OutputFlag = 0;
    grbparam.LogtoConsole = 0;
    sol = gurobi(model, grbparam);
    gamma = gamma * alpha_0;
    % Update iterate
    z_after = z_before + sol.x(batchtemp + 1:end);
    
    if norm(z_before - z_after) <= 1e-04
        break;
    end % End if
    
    % fprintf("Residual: %3.e \n", norm(x_before - x_after));
end % End for

x = x_after;
end % End function