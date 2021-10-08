function [sol, info] = proxsgdblindbatch(U, V, b, gamma, beta, init_z, maxiter, tol, ...
    early_stop, batch, adp_param, alpha_0, show_info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stop explaining and directly write code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get problem size
[m, n] = size(U);
adpmomentum = 0.01;

xidx = 1:n;
yidx = (n + 1) : (2 * n);

% Get initial point
z_before = init_z;
z_after = z_before;
bestz = init_z;

% Initialize arrays for storing results
% Array objs for maintaining the objective values
objs = zeros(maxiter * floor(m / batch) + 1, 1);
% Array best_onks for maintaining best objective values
bestobjs = zeros(maxiter * floor(m / batch) + 1, 1);

% Initialize trace related values
obj = sum(abs((U * z_before(xidx)) .* (V * z_before(yidx)) - b)) / m;
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
        nbatchiter = k * niter + idx;
        info.status = "Optimal";
        if early_stop
            if show_info
                disp("Optimizaition ends prematurely due to optimality");
            end % End if
            break;
        end
    end % End if
    
    idx = 0;
    
    for i = randperm(niter) % for i = randsample(1:m, m, true)
        idx = idx + 1;
        
        % pmt = randperm(m);
        % Sample from dataset
        batchidx = ((i - 1) * batch + 1) : (i * batch + resibatch * (i == niter));
        
        % Sample from dataset
        u = U(batchidx, :);
        v = V(batchidx, :);
        
        % Update momentum
        if beta == 999
            beta = adpmomentum / sqrt(k * m + idx);
            w = (1 + beta) * z_after - beta * z_before;
            beta = 999;
        else
            w = (1 + beta) * z_after - beta * z_before;
        end % End if
        
        % uTx and vTy are of size (batch * 1)
        uTx = u * z_after(xidx);
        vTy = v * z_after(yidx);
        
        z_before = z_after;
        
        sgn = sign(uTx .* vTy - b(batchidx));
        sgn = sgn + (2 * rand(batch + resibatch * (i == niter), 1) - 1).*(1 - abs(sgn));
        
        subgrad = mean([vTy .* u, uTx .* v] .* sgn, 1)';
        
        z_after = w - (alpha_0 * subgrad / gamma);
        
        obj = sum(abs((U * z_after(xidx)) .* (V * z_after(yidx)) - b)) / m;
        
        if obj < bestobj
            bestobj = obj;
            bestz = z_after;
        end % End if
        
        bestobjs((k * m - m) + idx + 1) = bestobj;
        objs((k * m - m) + idx + 1) = obj;
        
        if isnan(obj) || isinf(obj)
            info.status = "Diverged";
            break;
        end % End if
        
        if adp_param > 0
            % gamma = gamma + adp_param;
            gamma = sqrt(gamma^2 + adp_param);
        end % End if
        
    end % End for
    
    if isnan(obj) || isinf(obj)
        break;
    end % End if
    
    log = "- Epoch " + k + " - Obj: " + obj + " - Best obj: " + bestobj + ...
        " - Status: " + info.status;
    
    if show_info && mod(k, 50) == 0
        disp(log);
    end % End if
    
end % End for

% Collect information

% Solution array
sol.x = z_after;
sol.bestz = bestz;

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
end % End if

end % End function