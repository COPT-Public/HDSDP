function [x, y, s, kappa, tau] = lppotential(A, b, c, maxiter)
% Solve LP using potential reduction algorithm
% Minimize      c' * x
% Subject to     A * x == b
%                    x >= 0

[m, n] = size(A);

% Initialize
dim = 2 * m + 2 * n + 2;
rho = dim + sqrt(dim);

ypidx = 1:m; ymidx = m + 1 : 2 * m;
xidx = 2*m + 1 : 2*m + n; sidx = 2*m + h + 1 : 2*m +2*n;
kappaidx = dim - 1; tauidx = dim;
u_prev = ones(dim, 1) / dim;

% Complete the initial step using potential reduction
y = u(ypidx) - u(ymidx);
[f, g, pres, dres, compl] = getlpobjgrad(A, b, c, u_prev(xidx),...
                                         y, u_prev(sidx), u_prev(kappaidx), u_prev(tauidx));

potnew = rho * log(f) - sum(log(u_prev(xidx))) -...
                sum(log(u_prev(sidx))) -...
                log(u_prev(kappaidx)) -...
                log(u_prev(tauidx));
            
gk = rho / f * g - e./u_prev;

[potnew, gk] = potfg(rho, f, g, m, u_prev(xidx),...
                     u_prev(sidx), u_prev(kappaidx),...
                     u_prev(tauidx));

beta = 0.5;
lx = sum( (u_prev.^2).* gk * f ) / (norm(u_prev)^2 * rho);
px = u_prev.*(rho / f * (g - e * lx)) - e;
d = - beta / norm(px) * (u_prev .* px);
u_pres = u_prev + d;
potold = potnew;
recompute = true;


for i = 1:maxiter
    % Start potential reduction
    
    if recompute
        
        y = u(ypidx) - u(ymidx);
        [f, g, pres, dres, compl] = getlpobjgrad(A, b, c, u_pres(xidx),...
                                                y, u_pres(sidx), u_pres(kappaidx), u_pres(tauidx));
                                            
        mk = u_pres - u_prev;
        mk = mk / norm(mk);
        
        [potnew, gk] = potfg(rho, f, g, m, u_prev(xidx),...
                             u_prev(sidx), u_prev(kappaidx),...
                             u_prev(tauidx));
        
        % Gradient projection
        gk(m + 1:end) = gk(m + 1:end) - sum(gk(m + 1:end)) / (2 * n + 2);
        mk(m + 1:end) = mk(m + 1:end) - sum(mk(m + 1:end)) / (2 * n + 2);
        
        
        
        
    end % End if
    
    
    
    
    
    
    
    
end % End for


x = u_pres(xidx);
y = u_pres(ypidx) - u_pres(ymidx);
s = u_pres(sidx);
kappa = u_pres(kappaidx);
tau = u_pres(tauidx);

end % End function