function [lambda, delta, k] = dsdplanczosinv(Linv, dS)
% Compute the maximum eigenvalue of L^-1 dS * L^-T
if norm(dS, 'fro') < 1e-15
    lambda = 0;
    delta = 0;
end % End if

LinvT = Linv';
maxiter = 10;
% Initialization
[n, ~] = size(dS);
state = randn('state');
randn('state', 0);
v = randn(n, 1);
randn('state', state);

V = ones(n, maxiter + 1);
H = zeros(maxiter + 1, maxiter);
v = v / norm(v);
v = ones(n, 1);
v = v / norm(v);
V(:, 1) = v;

for k = 1:maxiter
    w = dS * (LinvT * v);
    w = Linv * w;
    
    wold = w;
    
    if (k > 1)
        w = w - H(k, k - 1) * V(:, k - 1);
    end % End if
    
    alp = w' * V(:, k);
    w = w - alp * V(:, k);
    H(k, k) = alp;
    
%     if norm(w) < 0.8 * norm(wold)
%         s = (w' * V(:, 1:k))';
%         w = w - V(:, 1:k) * s;
%         H(1:k, k) = H(1:k, k) + s;
%     end % End if
    
    nrm = norm(w);
    v = w / nrm;
    V(:, k + 1) = v;
    H(k + 1, k) = nrm;
    H(k, k + 1) = nrm;
    
    % Check convergence
    if mod(k, 5) == 0
        Hk = H(1:k, 1:k);
        Hk = 0.5 * (Hk + Hk');
        [Y, D] = eig(Hk);
        % Y = - Y;
        eigH = real(diag(D));
        [~, idx] = sort(eigH);
        res = abs(H(k+1, k) * Y(k, idx(k)));
        
        if (res <= 1e-04 || k == maxiter)
            lambda  = eigH(idx(k));
            lambda2 = eigH(idx(k - 1));
            z = V(:, 1:k) * Y(:, idx(k));
            z2 = V(:, 1:k) * Y(:, idx(k - 1));
            
            tmp  = dS * (LinvT * z);
            res1 = norm(Linv * tmp - lambda * z);
            tmp  = dS * (LinvT * z2);
            res2 = norm((Linv *  tmp) - lambda * z2);
            
            tmp = lambda - lambda2 - res2;
            if (tmp > 0)
                beta = tmp;
            else
                beta = eps;
            end % End if
            
            delta = min(res1, res1^2 / beta);
            if (delta <= 1e-03 || lambda + delta <= 0.8 )
                break;
            end % End if
        end % End if
        
    end % End if
end % End for

end % End function