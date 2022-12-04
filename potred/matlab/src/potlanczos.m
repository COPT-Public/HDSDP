function [zfin, lam, delta] = potlanczos(x, coneidx, rho, g, f, ATA, AT, A, scale, vstart)
% Compute minimum eigen-value using Lanczos iteration for projected Hessian
% matrix
%
%       H = P * [- (g * g') / f + ATA) + diag(d * f / rho)] * P
%
% The matrix-vector multiplication is specially processed

if nargin < 10
    vstart = [];
end % End if

rng(24);
n = size(AT, 1);
maxiter = 1000;

if isempty(vstart)
if scale
%     v = x;
%     v(coneidx) = 1.0;
    v = g;
    v(coneidx) = v(coneidx) .* x(coneidx);
else
    v = g;
end % End if
else
    v = vstart;
end % End if

% v = randn(n, 1);
V = zeros(n, maxiter + 1);
H = zeros(maxiter + 1, maxiter);
v = v / norm(v);
V(:, 1) = v;
tol = 1e-06;
ncone = length(coneidx);
useold = false;

for k = 1:maxiter
    
    %     fprintf("%f \n", v' * Hproj * v);
    % w = P * z
    if useold
        w = v;
        w(coneidx) = w(coneidx) - sum(w(coneidx)) / ncone;
        w1 = (- (g' * w) / f) * g;
        w2 = AT * (A * w); % Or w2 = ATA * w;
        w3 = (f / rho) * (d .* w);
        w = - w1 - w2 - w3;
        w(coneidx) = w(coneidx) - sum(w(coneidx)) / ncone;
    else
        w = -MXv(v, coneidx, AT, A, ATA, x, f, g, rho, scale);
    end % End if
    
    wold = w;
    
    if (k > 1)
        w = w - H(k, k-1) * V(:, k-1);
    end % End if
    
    alp = w' * V(:, k);
    w = w - alp * V(:, k);
    H(k, k) = alp;
    
    if (norm(w) <= 0.99 * norm(wold) || 1)
        s = (w' * V(:,1:k))';
        w = w - V(:,1:k) * s;
        H(1:k, k) = H(1:k, k) + s;
    end % End if
    
    nrm = norm(w);
    v = w / nrm;
    V(:, k+1) = v;
    H(k+1,k) = nrm;  H(k,k+1) = nrm;
    
    if (rem(k, 5) == 0 || k == maxiter)
        
        Hk = H(1:k,1:k); Hk = 0.5 * (Hk + Hk');
        [Y, D] = eig(Hk);
        eigH  = real(diag(D));
        [~,idx] = sort(eigH);
        res_est = abs(H(k+1, k) * Y(k, idx(k)));
        
        if (res_est <= 0.01) || (k == maxiter)
            lam = eigH(idx(k));
            lam2 = eigH(idx(k-1));
            z = V(:, 1:k) * Y(:, idx(k));
            z2 = V(:, 1:k) * Y(:, idx(k - 1));
            
            if useold
                zz = z;
                zz(coneidx) = zz(coneidx) - sum(zz(coneidx)) / ncone;
                zz1 = (- (g' * zz) / f) * g;
                zz2 = AT * (A * zz);
                zz3 = (f / rho) * (d .* zz);
                zz = - zz1 - zz2 - zz3;
                zz(coneidx) = zz(coneidx) - sum(zz(coneidx)) / ncone;
            else
                zz = -MXv(z, coneidx, AT, A, ATA, x, f, g, rho, scale);
            end % End if
            
            res = norm(zz - lam * z);
            zfin = zz;
            
            if useold
                zz = z2;
                zz(coneidx) = zz(coneidx) - sum(zz(coneidx)) / ncone;
                zz1 = (- (g' * zz) / f) * g;
                zz2 = AT * (A * zz);
                zz3 = (f / rho) * (d .* zz);
                zz = - zz1 - zz2 - zz3;
                zz(coneidx) = zz(coneidx) - sum(zz(coneidx)) / ncone;
            else
                zz = -MXv(z2, coneidx, AT, A, ATA, x, f, g, rho, scale);
            end % End if
            
            res2 = norm(zz - lam * z2);
            tmp = lam - lam2 - res2;
            
            if tmp > 0
                beta = tmp;
            else
                beta = eps;
            end % End if
            
            delta = min(res, res^2 / beta);
            
            if delta <= tol || (lam > 0 && lam - abs(delta) > 0)
                fprintf("Lanzos iteration %d \n", k);
                break;
            end % End if
            
        end % End if
    end % End if
end % End if

end % End function