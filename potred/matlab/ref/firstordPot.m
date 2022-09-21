function [x] = firstordPot(A, gamma)
% Implement the potential reduction algorithm for simplex-constrained QP

[~, n] = size(A);

rho = n + sqrt(n);
e = ones(n, 1); 
x = e / n;
ATA = A' * A;

potold = inf;
potbeta = 0.5;

for i = 1:10000
    
%     fprintf("%10.3e %10.3e \n", x(1), x(2));
    
    [f, g] = fpot(A, ATA, x);
    potnew = rho * log(f) - sum(log(x));
    potred = potnew - potold;
    gk = rho / f * g - e./x;
    
%     beta = max(0.1, 1 / (2 + rho * gamma / f));
    if potred < 0
        beta = potbeta; % 1 / (2 + rho * gamma / f);
    else
        potbeta = potbeta / 2;
        beta = 1 / (2 + rho * gamma / f);
    end % End if
    
    lx = sum( (x.^2).* gk * f ) / (norm(x)^2 * rho);
    px = x.*(rho / f * (g - e * lx)) - e;
    d = - beta / norm(px) * (x .* px);
    x = x + d;
    
    potold = potnew;
%     fprintf("%6.3e  %6.3e\n", f, potnew);
    fprintf("%3d  %6.3e  %6.3e\n", i, f, potred);
    
end % End for

assert(min(x) >= 0);


end % End function