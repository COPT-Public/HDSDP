function [sol, info] = proxptblind(U, V, b, gamma, beta, init_z, maxiter, tol, ...
    early_stop, adp_param, alpha_0, show_info)

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
objs = zeros(m * maxiter + 1, 1);
% Array best_onks for maintaining best objective values
bestobjs = zeros(m * maxiter + 1, 1);

% Initialize trace related values
obj = sum(abs((U * z_before(xidx)) .* (V * z_before(yidx)) - b)) / m;
bestobj = obj;
objs = objs + obj;
bestobjs= bestobjs + bestobj;

% Number of epochs before reaching tolerance
nepochs = maxiter;
nbatchiter = maxiter * m;

% Initialize info struct
info.status = "Not Optimal";

for k = 1:maxiter
    
    if bestobj < tol && nepochs == maxiter
        nepochs = k;
        nbatchiter = k * m + idx;
        info.status = "Optimal";
        if early_stop
            if show_info
                disp("Optimizaition ends prematurely due to optimality");
            end % End if
            break;
        end
    end % End if
    
    idx = 0;
    niter = m;
    
    for i = randperm(niter) % for i = randsample(1:m, m, true)
        idx = idx + 1;
        
        gamma = gamma / alpha_0;
        
        % Sample from dataset
        u = U(i, :);
        v = V(i, :);
                
        % Update momentum
        if beta == 999
            beta = 1 / alpha_0 * adpmomentum / sqrt(k * m + idx);
            w = (1 + beta) * z_after - beta * z_before;
            beta = 999;
        else
            w = (1 + beta) * z_after - beta * z_before;
        end % End if
        
        uTx = u * z_after(xidx);
        vTy = v * z_after(yidx);
        
        z_before = z_after;
        
        % Do proximal point update
        vTwy = v * w(yidx);
        uTwx = u * w(xidx);
        unorm = norm(u);
        vnorm = norm(v);
        
        if abs(uTx * vTy - b) >= 1e-08
            
            p1 = w;
            p1(xidx) = p1(xidx) - (gamma * vTwy - vnorm^2 * uTwx) / (gamma^2 - norm(u)^2 * norm(v)^2) * u';
            p1(yidx) = p1(yidx) - (gamma * uTwx - unorm^2 * vTwy) / (gamma^2 - norm(u)^2 * norm(v)^2) * v';
            
            p2 = w;
            p2(xidx) = p2(xidx) - (- gamma * vTwy - vnorm^2 * uTwx) / (gamma^2 - norm(u)^2 * norm(v)^2) * u';
            p2(yidx) = p2(yidx) - (- gamma * uTwx - unorm^2 * vTwy) / (gamma^2 - norm(u)^2 * norm(v)^2) * v';
            
            p3 = w;
            p3(xidx) = p3(xidx) - (gamma * vTwy - vnorm^2 * uTwx) / (gamma^2 - norm(u)^2 * norm(v)^2) * u';
            p3(yidx) = p3(yidx) - (- gamma * uTwx - unorm^2 * vTwy) / (gamma^2 - norm(u)^2 * norm(v)^2) * v';
            
            p4 = w;
            p4(xidx) = p4(xidx) - (- gamma * vTwy - vnorm^2 * uTwx) / (gamma^2 - norm(u)^2 * norm(v)^2) * u';
            p4(yidx) = p4(yidx) - (gamma * uTwx - unorm^2 * vTwy) / (gamma^2 - norm(u)^2 * norm(v)^2) * v';
            
            fvals = [abs((u * p1(xidx)) * (v * p1(yidx)) - b(i)) + (gamma / 2) * norm(p1 - w)^2, ...
                     abs((u * p2(xidx)) * (v * p2(yidx)) - b(i)) + (gamma / 2) * norm(p2 - w)^2, ...
                     abs((u * p3(xidx)) * (v * p3(yidx)) - b(i)) + (gamma / 2) * norm(p3 - w)^2, ...
                     abs((u * p4(xidx)) * (v * p4(yidx)) - b(i)) + (gamma / 2) * norm(p4 - w)^2];
            [~, minidx] = min(fvals);
            
            if minidx == 1
                z_after = p1;
            elseif minidx == 2
                z_after = p2;
            elseif minidx == 3
                z_after = p3;
            else 
                z_after = p4;
            end% End if
            
        else
            rts = roots([vnorm^2, - vnorm^2 * uTwx, 0, b(i) * unorm^2 * vTwy, - b(i)^2 * unorm^2]);
            fval = inf;
            for q = 1:4
                eta = rts(q);
                zeta = (eta * uTwx - eta^2) / (b(i) * unorm^2);
                
                pt = w;
                
                pt(xidx) = pt(xidx) - zeta * b(i) / eta * u';
                pt(yidx) = pt(yidx) - zeta * eta * v';
                
                if fval > abs((u * pt(xidx)) * (v * pt(yidx)) - b(i)) + gamma / 2 * norm(pt - w)^2
                    z_after = pt;
                end % End if
                    
            end % End for
            
        end % End if
        
        obj = sum(abs((U * z_after(xidx)) .* (V * z_after(yidx)) - b)) / m;
        
        gamma = gamma * alpha_0;
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