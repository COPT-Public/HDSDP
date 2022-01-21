function [Ry, S, y] = dsdpInitialize(A, C, tau, dsdpInitializeStrategy, initbeta)
% Initializer for DSDP6. Get Ry such that
% S = C * tau - Ry >= 0
% We note that y is always initialized by 0

if nargin < 3
    tau = 1;
end % End if

[n, ~] = size(C);
m = length(A);

y = zeros(m, 1);
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
    S = Ctau - Ry;
    
    while ~ dsdpIspsd(S)
        alpha = alpha * 2;
        S = Ctau - Ry;
    end % End while
    
    S = S + initbeta * eyemat;
    
elseif dsdpInitializeStrategy == "fro"
    
    if norm(Ctau, 'fro') == 0
        Ry = - speye(n) * initbeta;
    else
        Ry = - speye(n) * (norm(Ctau, 'fro')) * initbeta;
    end % End if 
    
    S  = Ctau - Ry; 
else
    S  = speye(n);
    Ry = Ctau - S;
end % End if
    
end % End function