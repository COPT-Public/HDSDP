function [D, E, Ascal] = ruizscale2(A, maxiter)
% Implement the Ruiz scaling algorithm for general matrix

if nargin == 1
    maxiter = 100;
end % End if

[m, n] = size(A);

B = A;

d1 = ones(m, 1);
d2 = ones(n, 1);

for i = 1:maxiter
    
    BB = B.^2;
    d1 = d1 ./ (sqrt(sum(BB, 2)).^0.5);
    d2 = (m / n)^0.25 * (d2 ./ (sqrt(sum(BB, 1)').^0.5));
    
    D1 = diag(d1);
    D2 = diag(d2);
    
    B = D1 * A * D2;
    
    if min(d1) <= 1e-06 * max(d1) && min(d2) <= 1e-06 * max(d2)
        break;
    end % End if
    
end % End for

D = d1;
E = d2;
Ascal = A;

end % End function