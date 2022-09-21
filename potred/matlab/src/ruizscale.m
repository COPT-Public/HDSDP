function [D, E, Ascal] = ruizscale(A, maxiter)
% Implement the Ruiz scaling algorithm for general matrix

if nargin == 1
    maxiter = 100;
end % End if

[m, n] = size(A);
D = ones(m, 1);
E = ones(n, 1);

for i = 1:maxiter
    
    dR = sqrt(max(abs(A), [], 2));
    dC = max(abs(A)).^(-1/2);
    R = diag(dR);
    C = diag(dC);
    A = R \  (A * C);
    D = D ./ dR;
    E = E .* dC';
    
    if norm(dR - 1, 'inf') <= 1e-08 && norm(dC - 1, 'inf') <= 1e-08
        break;
    end % End if
    
end % End for

Ascal = A;

end % End function