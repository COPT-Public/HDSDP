function [sol, info] = proxptbatch(A, b, gamma, init_x, maxiter, tol, ...
    early_stop, batch, adp_param, alpha_0, show_info, beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Inertial Model-based Stochastic Methods                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Oct 7th, 2020                                %
%                      Phase Retrieval Problem                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements the optimization algorithm with respect to    %
% proximal point method with momentum and variable metric to solve       %
% Phase Retrieval Problem                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                                 %
%          A, b : data in optimization                                   %
%   gamma, beta : optimization parameters, gamma for stepsize parameter  %
%                 and beta is for momentum parameter                     %
%        init_x : initial starting point of algorithm                    %
%       maxiter : maximum number of iterations (epochs) allowed          %
%           tol : tolerance allowed to end the algorithm prematurely     %
%    early_stop : whether to stop optimization when tolerance is reached %
%     adp_param : parameters for adaptive method                         %
%        use_vm : parameters for variable metric                         %
%     show_info : whether to display optimization progress               %
% Output:                                                                %
%           sol : an array storing the last/best solution after the      %
%                 optimization procedure                                 %
%          info : an array storing the information related to the        %
%                 optimization procedure                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 11
    beta = 0;
end % End if 

% Get problem size
[m, n] = size(A);

% Get initial point
x_before = init_x;
x_after = x_before;
bestx = init_x;

% Array objs for maintaining the objective values
objs = zeros(maxiter * floor(m / batch) + 1, 1);
% Array best_onks for maintaining best objective values
bestobjs = zeros(maxiter * floor(m / batch) + 1, 1);

% Initialize trace related values
obj = sum(abs((A * x_before).^2 - b)) / m;
bestobj = obj;
objs = objs + obj;
bestobjs= bestobjs + bestobj;

% Number of epochs before reaching tolerance
nepochs = maxiter;
nbatchiter = maxiter * m / batch;

% Initialize info struct
info.status = "Not Optimal";
if mod(m, batch) == 0
        niter = m / batch;
        resibatch = 0;
    else
        niter = floor(m / batch);
        resibatch = mod(m, batch);
end % End if

for k = 1:maxiter
    
    if bestobj < tol && nepochs == maxiter
        nepochs = k;
        nbatchiter = k * m / batch + idx;
        info.status = "Optimal";
        if early_stop
            if show_info
                disp("Optimizaition ends prematurely due to optimality");
            end % End if
            break;
        end
    end % End if
    
    idx = 0;
    
    %     if adp_param > 0
    %         % gamma = gamma + adp_param;
    %         gamma = sqrt(gamma^2 + k * m * adp_param);
    %     end % End if
       
    for i = randperm(niter) % for i = randsample(1:m, m, true)
        idx = idx + 1;
        
        % Sample from dataset
       % Sample from dataset
        batchidx = ((i - 1) * batch + 1) : (i * batch + resibatch * (i == niter));
        a = A(batchidx, :);
        
        batchtemp = batch + resibatch * (i == niter);
        y = (1 + beta) * x_after - beta * x_before;
        
        % Update momentum
        x_before = x_after;        
        gamma = gamma / alpha_0;
        
        x_after = proxlinsolve(a, b(batchidx), y, gamma, 1, 50);
        
%         % Solve mini-batch sub-problem using QCQP
%         qpn = batch + n;
%        
%         % Construct objective coefficient
%         model.obj = zeros(qpn, 1);
%         model.obj(1:batch) = 1 / batch;
%        
%         % Get lowerbound
%         model.lb = - inf(qpn, 1);
%         
%         % Get A and rhs
%         model.A = sparse(1, qpn);
%         
%         % Get quadratic constraints
%         for j = 1:batch
%             
%             % Get quadratic matrix
%             Q0 = sparse(batch, batch);
%             normx = gamma * norm(x_after)^2 / 2;
%             Q1 = speye(n) * gamma / 2 + a(j, :)' * a(j, :);
%             Q2 = speye(n) * gamma / 2 - a(j, :)' * a(j, :);
%             model.quadcon(2 * j - 1).Qc = [Q0, sparse(batch, n);
%                                            sparse(n, batch), Q1];
%             model.quadcon(2 * j - 1).q = [zeros(batch, 1); - gamma * x_before];
%             model.quadcon(2 * j - 1).q(j) = -1;
%             model.quadcon(2 * j - 1).rhs = b(batchidx(j)) - normx;
%             
%             model.quadcon(2 * j).Qc = [Q0, sparse(batch, n);
%                                            sparse(n, batch), Q2];
%             model.quadcon(2 * j).q = [zeros(batch, 1); - gamma * x_before];
%             model.quadcon(2 * j).q(j) = -1;
%             model.quadcon(2 * j).rhs = - b(batchidx(j)) - normx;
%             
%         end % End for
%         
%         % Solve model
%         grbparam.OutputFlag = 0;
%         grbparam.LogtoConsole = 0;
%         sol = gurobi(model, grbparam);
%         x_after = sol.x(batch + 1 : end);
        
        gamma = gamma * alpha_0;
        
        obj = sum(abs((A * x_after).^2 - b)) / m;
        
        if obj < bestobj
            bestobj = obj;
            bestx = x_after;
        end % End if
        
        bestobjs((k * niter - niter) + idx + 1) = bestobj;
        objs((k * niter - niter) + idx + 1) = obj;
        
        if isnan(obj)
            info.status = "Diverged";
        end % End if
        
        if adp_param > 0
            % gamma = gamma + adp_param;
            gamma = sqrt(gamma^2 + adp_param);
        end % End if
        
    end % End for
    
    log = "- Epoch " + k + " - Obj: " + obj + " - Best obj: " + bestobj + ...
        " - Status: " + info.status;
    
    if show_info && mod(k, 1) == 0
        disp(log);
    end % End if
    
end % End for

% Collect information
% Solution array
sol.x = x_after;
sol.bestx = bestx;

% Information array
info.nepochs = nepochs;
info.niter = nbatchiter;
info.objs = objs;
info.bestobjs = bestobjs;

% Display summary
if show_info && info.status == "Optimal"
    disp("- Algorithm reaches optimal after " + nepochs + " epochs (" + ...
        nbatchiter + " iterations)");
elseif show_info && info.status == "Not Optimal"
    disp("- Algorithm fails to reach desired accuracy after " +...
        nepochs + " epochs");
elseif show_info && info.status == "Diverged"
    disp("- Algorithm diverges");

end % End function

