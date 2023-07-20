clear;

printconestat = 1;
printdualslack = 0;
setupkkt = 1;

[At, b, c, K] = readsdpa(fullfile('/Users/gaowenzhi/Desktop/gwz/benchmark/sdplib', 'mcp100.dat-s'));

try
    [nsqr, ~] = size(At);
    A = At';
catch
    [~, nsqr] = size(A);
end % End try

l = K.l;
s = length(K.s);
[m, ~] = size(A);
Amat = cell(m, s);
Cmat = cell(1, s);

counter = 1 + l;
for q = 1:s
    n = K.s(q);
    Cmat{q} = reshape(c(counter:counter + n * n - 1), n, n);
    counter = counter + n * n;
end % End for

Alp = A(:, 1:l);
for i = 1:m
    counter = 1 + l;
    for q = 1:s
        n = K.s(q);
        Amat{i, q} = reshape(A(i, counter:counter + n * n - 1), n, n);
        counter = counter + n * n;
    end % End for
end % End for

% Count block statistics
for q = 1:s
    
    n = K.s(q);
    zero = 0;
    sps = 0;
    ds = 0;
    
    if printconestat
        fprintf("| Block %d \n", q);
        fprintf("| %4d | %10.3e %10.3e \n", 0, absnorm(Cmat{1, q}), norm(Cmat{1, q}, 'fro'));
    end % End if
    
    for i = 1:m
        Acoef = Amat{i, q};
        if printconestat && absnorm(Acoef) > 0.0
            fprintf("| %4d | %10.3e %10.3e \n", i, absnorm(Acoef), norm(Acoef, 'fro'));
        end % End if
        if length(nonzeros(tril(Acoef))) == 0 %#ok
            zero = zero + 1;
        elseif length(nonzeros(tril(Acoef))) < 0.3 * (n * (n + 1)) / 2
            sps = sps + 1;
        else
            ds = ds + 1;
        end % End if
    end % End for
    
    if printconestat
        fprintf("Zero: %d Nonzero: %d Sparse: %d Dense: %d\n",...
            zero, m - zero, sps, ds);
    end % End if
    
    if printdualslack
        Rd = 1e+06;
        y = 1:m;
        B = eye(n) * Rd - hdsdp_aty(Amat(:, q), y) + Cmat{q} * 1.5;
        spB = sparse(B);
        fprintf("logdet: %20.10e \n", hdsdp_logdet(spB));
    end % End if
    
    if setupkkt
        Rd = 1000.0;
        y = 1:m;
        B = eye(n) * Rd - hdsdp_aty(Amat(:, q), y) + Cmat{q} * 1.5;
        [M, asinv, asinvrdsinv, asinvcsinv, csinv, csinvcsinv, csinvrdcsinv] =...
        hdsdp_kktbuild(Amat, Cmat{1}, B, Rd);
    end % End if
    
end % End for

