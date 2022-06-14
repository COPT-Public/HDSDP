function [Ry, S, y] = dsdpInitialize(A, C, tau, dsdpInitializeStrategy, initbeta, y0)
% Initializer for DSDP6. Get Ry such that
% S = C * tau - Ry >= 0
% We note that y is always initialized by 0

if nargin < 3
    tau = 1;
end % End if

[n, ~] = size(C);
m = length(A);

y = zeros(m, 1);

if ~isempty(y0)
    y = y0;
end % End if
    
eyemat = speye(n, n);

if tau ~= 1
    Ctau = C * tau;
else
    Ctau = C;
end % End if

if length(unique(C)) == 1
    initbeta = initbeta * 10;
end % End if

% A choice using minimum eigenvalue
if dsdpInitializeStrategy == "eigs"
    
    Ry = 10 * eyemat * eigs(Ctau, 1, 'smallestreal');
    S = Ctau - Ry;
    
elseif dsdpInitializeStrategy == "linesearch"
    
    alpha = 1.0;
    Ry = - alpha * eyemat;
        
    if isempty(y0)
        S = Ctau - Ry;
        
        while ~ dsdpIspsd(S)
            alpha = alpha * 2;
            Ry = alpha * Ry;
            S = Ctau - Ry;
        end % End while
        S = S + 0.0 * eyemat;
    else
       S = Ctau - dsdpgetATy(A, y0) - Ry;
       while ~ dsdpIspsd(S)
            alpha = alpha * 2;
            S = Ctau - Ry;
            Ry = alpha * Ry;
        end % End while
    end % End if
    
elseif dsdpInitializeStrategy == "fro"
    
    if norm(Ctau, 'fro') == 0
        Ry = - speye(n) * initbeta;
    else
        Ry = - speye(n) * max(norm(Ctau, 'fro'), 100) * initbeta;
    end % End if 
    
    S  = Ctau - Ry; 
else
    S  = speye(n);
    Ry = Ctau - S;
end % End if
    
end % End function