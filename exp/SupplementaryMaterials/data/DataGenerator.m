kappa_range = [1, 10];
pfail_range = [0, 0.2, 0.3];

optx = randn(100, 1);
optx = optx / norm(optx);


for kappa = kappa_range
    for pfail = pfail_range
        clear data;
        [A, b] = GetData(optx, 300, kappa, pfail);
        b = max(b, 0);
        data.A = A;
        data.b = b;
        data.optx = optx;
        data.bestloss = sum(abs((A * optx).^2 - b)) / 300;
        save("kappa_" + kappa + "_pfail_" + pfail * 10 + ".mat", "data");
        
        clear data;
        [U, V, b] = GetDataBlind(optx, 300, kappa, pfail);
        data.U = U;
        data.V = V;
        data.b = b;
        data.optx = optx;
        data.bestloss = sum(abs((U * optx) .* (V * optx) - b)) / 300;
        save("blind_kappa_" + kappa + "_pfail_" + pfail * 10 + ".mat", "data");
        
    end % End for
end % End for

% % Generate initial x
% init_x = randn(100, 1);
% save("init_x.mat", "init_x");