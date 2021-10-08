function [sol, info] = proxlin(A, b, gamma, beta, init_x, maxiter, tol, ...
    early_stop, adp_param, alpha_0, use_vm, show_info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Inertial Model-based Stochastic Methods.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Oct 7th, 2020                                %
%                      Phase Retrieval Problem                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements the optimization algorithm with respect to    %
% proximal linear method with momentum and variable metric to solve      %
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

% Get problem size
[m, n] = size(A);
adpmomentum = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create backdoor for images
imgcount = 999;
imgres = 0;
if m == 768
    imgres = zeros(11, 256);
    imgres(1, :) = init_x';
    imgcount = 1;
end % End if
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get initial point
x_before = init_x;
x_after = x_before;
bestx = init_x;

% Initialize arrays for storing results
% Array objs for maintaining the objective values
objs = zeros(m * maxiter + 1, 1);
% Array best_onks for maintaining best objective values
bestobjs = zeros(m * maxiter + 1, 1);

% Initialize trace related values
obj = sum(abs((A * x_before).^2 - b)) / m;
bestobj = obj;
objs(1) = obj;
bestobjs(1) = bestobj;

% Number of epochs before reaching tolerance
nepochs = maxiter;
niter = maxiter * m;

% Variable metric vector, starting from identity
metric = ones(n, 1);
if use_vm
    metric = metric * 0.1;
end % End

% Initialize info struct
info.status = "Not Optimal";

for k = 1:maxiter
    
    if bestobj < tol && nepochs == maxiter
        nepochs = k;
        niter = k * m + idx;
        info.status = "Optimal";
        if early_stop
            if show_info
                disp("Optimizaition ends prematurely due to optimality");
            end % End if
            break;
        end
    end % End if
    
    idx = 0;
    vm_update = false;
    
    %     if adp_param > 0
    %         % gamma = gamma + adp_param;
    %         gamma = sqrt(gamma^2 + k * m * adp_param);
    %     end % End if
    
    for i = randperm(m) % for i = randsample(1:m, m, true)
        idx = idx + 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mod(k * m + idx, 59) == 0 && imgcount <= 10
            imgcount = imgcount + 1;
            imgres(imgcount, :) = bestx;
        end % End if 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Sample from dataset
        a = A(i, :);
        
        % Update momentum
        if beta == 999
            beta = 1 / alpha_0 * adpmomentum / sqrt(k * m + idx);
            y = (1 + beta) * x_after - beta * x_before;
            beta = 999;
        else
            y = (1 + beta) * x_after - beta * x_before;
        end % End if 
           
        gamma = gamma / alpha_0;
        x_before = x_after;
        aTx = a * x_before;
        ksi = aTx^2 + 2 * aTx * a * (y - x_after) - b(i);
        zeta = 2 * aTx * a';
        coef = - ksi / ((zeta./metric)'*zeta) * gamma;
        coef = sign(coef) * min(1, abs(coef));
        x_after = y + (coef * zeta / gamma)./metric;
        gamma = gamma * alpha_0;
        obj = sum(abs((A * x_after).^2 - b)) / m;
        
        if obj < bestobj
            bestobj = obj;
            bestx = x_after;
        end % End if
        
        bestobjs(k * m - m + idx + 1) = bestobj;
        objs(k * m - m + idx + 1) = obj;
        
        if isnan(obj)
            info.status = "Diverged";
        end % End if
        
        % Update metric across iteration
        if use_vm
            % Compute metric incremental
            sgn = sign(aTx^2 - b(idx));
            
            if sgn
                subgrad = 2 * aTx * a' * sgn;
            else
                subgrad = 2 * aTx * a' * (2 * rand() - 1);
            end % End if
            
            vm_increment = max((sqrt(diag(subgrad * subgrad')) - metric), 0) /...
                (idx + k * m - m + 1);
            
            metric = sqrt(metric.^2 + use_vm * diag(subgrad * subgrad'));
            % metric = metric + vm_increment * use_vm;
            % metric = sqrt(metric.^2 + vm_increment * use_vm);
            
            if sum(vm_increment) > 0
                vm_update = true;
            end % End if
            
        end % End if
        
        if adp_param > 0
            % gamma = gamma + adp_param;
            gamma = sqrt(gamma^2 + adp_param);
        end % End if
        
    end % End for
    
    log = "- Epoch " + k + " - Obj: " + obj + " - Best obj: " + bestobj + ...
        " - Status: " + info.status;
    
    if use_vm && vm_update
        log = log + " (metric updated)";
    end % End if
    
    if show_info && mod(k, 50) == 0
        disp(log);
    end % End if
    
end % End for
% Collect information

% Solution array
sol.x = x_after;
sol.bestx = bestx;

% Information array
info.nepochs = nepochs;
info.niter = niter;
info.objs = objs;
info.bestobjs = bestobjs;
info.imgres = imgres;

% Display summary
if show_info && info.status == "Optimal"
    disp("- Algorithm reaches optimal after " + nepochs + " epochs (" + ...
        niter + " iterations)");
elseif show_info && info.status == "Not Optimal"
    disp("- Algorithm fails to reach desired accuracy after " +...
        nepochs + " epochs");
elseif show_info && info.status == "Diverged"
    disp("- Algorithm diverges");
end % End if

end % End function

