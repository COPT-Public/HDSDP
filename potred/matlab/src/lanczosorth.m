function [v, e] = lanczosorth(A)

maxiter = 30;
n = size(A, 1);
EvalFLSO = zeros(maxiter, maxiter);
ErrorBoundsFLSO = zeros(m,m);
RitzComponentsFLSO = zeros(m,m);
RitzComponentsFLSO = zeros(maxiter, maxiter);
SelectFLSO = 0;

q = rand(n, 1);
qq = q / norm(q);
QFLSO = zeros(n, maxiter + 1);
QFLSO(:, 1) = qq;

for i = 1:maxiter
    
    z = A * qq;
    alpha(i) = qq' * z;
    z = z - alpha(i) * qq;
    
    if i > 1
        z = z - beta(i - 1) * QFLSO(:, i - 1);
    end % End if
    
    beta(i) = norm(z);
    qq = z / beta(i);
    
    if i > 0
        Ti = diag(alpha(1:i)) + diag(beta(1:i-1),1) + diag(beta(1:i-1),-1);
    else
        Ti = alpha(1);
    end % End if
    
    [V,D] = eig(Ti);
    
    [Ds,Is] = sort(-diag(D));
    EvalFLSO(1:i,i) = -Ds;
    Vs = V(:,Is);
    ErrorBoundsFLSO(1:i,i) = abs((beta(i)*Vs(i,:))');
    select = find( ErrorBoundsFLSO(1:i,i) <= sqrt(eps)*max(abs(Ds)) );
    
    if (length(select) > 0)
        SelectFLSO(1:length(select),i) = select;
        for s = select'
            ritzvector = QFLSO*Vs(:,s);
            z = z - (ritzvector'*z)*ritzvector;
        end
        beta(i) = norm(z);
        qq = z/beta(i);
    else
        SelectFLSO(1,i) = 0;
    end
    
    QFLSO(:, i + 1) = qq;
    RitzComponentsFLSO(1:i,i) = ((QFLSO(:, i+1)' * QFLSO(:,1:i)) * Vs)';
    
end % End for

[qQ,rQ] = qr(QFLSO(:,1:m),0);

end % End function