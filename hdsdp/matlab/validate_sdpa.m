clear;

[At, b, c, K] = readsdpa(fullfile('/Users/gaowenzhi/Desktop/dsdp6/hsd/benchmark/sdplib', 'arch8.dat-s'));

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

for q = 1:s
    counter = 1 + l;
    n = K.s(q);
    Cmat{q} = reshape(c(counter:counter + n * n - 1), n, n);
    counter = counter + n * n;
end % End for

Alp = A(:, 1:l);
for i = 1:m
    counter = 1 + l;
    for q = 1:s
        n = K.s(q);
        Amat{i, q} = full(reshape(A(i, counter:counter + n * n - 1), n, n));
        counter = counter + n * n;
    end % End for
end % End for

% Count block statistics
for q = 1:s
    
    n = K.s(q);
    zero = 0;
    sps = 0;
    ds = 0;
    
    fprintf("| Block %d \n", q);
    fprintf("| %4d | %10.3e %10.3e \n", 0, absnorm(Cmat{1, q}), norm(Cmat{1, q}, 'fro'));
    for i = 1:m
        Acoef = Amat{i, q};
        fprintf("| %4d | %10.3e %10.3e \n", i, absnorm(Acoef), norm(Acoef, 'fro'));
        if length(nonzeros(tril(Acoef))) == 0 %#ok
            zero = zero + 1;
        elseif nonzeros(tril(Acoef)) < 0.3 * (n * (n + 1)) / 2
            sps = sps + 1;
        else
            ds = ds + 1;
        end % End if
    end % End for
    
    fprintf("Zero: %d Nonzero: %d Sparse: %d Dense: %d\n",...
        zero, m - zero, sps, ds);
    
end % End for
