clear;
clc;

rng(24);

warning off;
[At,b,c,K] = readsdpa(fullfile('.', 'sdplib', 'hamming_9_5_6.dat-s'));

try
    [nsqr, ~] = size(At);
    A = At';
catch
    [~, nsqr] = size(A);
end % End try

l = K.l;
s = length(K.s);
[m, ~] = size(A);
Amat = cell(m, 1);

for i = 1:m
    counter = 1;
    B = sparse(0, 0);
    n = K.l;
    B = blkdiag(B, diag(A(i, counter:counter + n - 1)));
    counter = counter + n;
    for q = 1:s
        n = K.s(q);
        B = blkdiag(B, reshape(A(i, counter:counter + n * n - 1), n, n));
        counter = counter + n * n;
    end % End for
    Amat{i} = B;
end % End for

counter = 1;
n = K.l;
B = sparse(0, 0);
B = blkdiag(B, diag(c(counter:counter + n - 1)));
counter = counter + n;
for q = 1:s
    n = K.s(q);
    B = blkdiag(B, reshape(c(counter:counter + n * n - 1), n, n));
    counter = counter + n * n;
end % End for

C = B;

rng(24);

if false
    y = zeros(m, 1);
    for i = 1:500
        [emax, grad] = eigrad(Amat, C, y);
        y = y - 15 * grad / sqrt(i);
        if mod(i, 10) == 0 || i == 1
            fprintf("%d  %e \n", i, emax);
        end % End if
        if emax <= 0.0
            fprintf("%d  %e \n", i, emax);
            break;
        end % End if
    end % End for
end % End if

tol = 1e-04;
y0 = zeros(m, 1);
[emax, grad] = eigrad(Amat, C, y0);
fprintf("%d %e \n", 0, emax);


total = 1e-04;
fbest = emax;
baserate = 1.1;
lambda = 0.1;
maxiter = max(20, min(m / 10, 100));
baseiter = max(maxiter / 10, 20);

for i = 1:1000
    
    [y0, fub, iter] = dsdpAPLsolve(y0, Amat, C, baseiter, baserate * 2 * i / (2 * i + 1), lambda);
    total = total + iter;
    
    if (fub + lambda) > (fbest + lambda) * 2 * i / (2 * i + 1)
        baseiter = ceil(baseiter * 1.1);
        baserate = min(0.9,  baserate * 1.1);
    elseif (fub + lambda) < (fbest + lambda) * 2 * i / (2 * i + 1) * 0.8
        baserate = baserate * 0.9;
        
        if lambda > 0
            lambda = lambda * 0.9;
        else
            lambda = lambda / 0.9;
        end % End if 
        
    end % End if
    
    if fub < tol
        lambda = 1;
    end % End if
    
    fbest = min(fub, fbest);
    
    if (fub + lambda) < tol || total >= maxiter
        break;
    end % End if
    
    fprintf("%2d %10.3e %10.3e\n", i, fub, -lambda);
    
end % End for

fprintf("%e \n", fub);



