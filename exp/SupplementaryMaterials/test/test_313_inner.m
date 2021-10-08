% Inner test script for Synthetic Phase Retrieval Momentum&Minibatch (3.1.3)
A = data.A;
b = data.b;
bestloss = data.bestloss;

tol = bestloss * 1.5;

[m, n] = size(A);

nSgdEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nSgdmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nSgdbEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nSgdbmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxLinEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxLinmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxLinbEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxLinbmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxPtEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxPtmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxPtbEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxPtbmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;

for k = 1:length(steprange)
    
    stepsize = steprange(k);
    
    parfor i = 1:nTest
        
        init_x = randn(n, 1);
        
        % Nothing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test SGD
        [sgdsol, sgdinfo] = proxsgd(A, b, sqrt(maxiter * m), 0, init_x, ...
            maxiter, tol, true, 1, 0, stepsize, 0, show_info);
        
        if sgdinfo.status == "Optimal"
            nSgdEpochtoOpt(i, k) = sgdinfo.niter;
        end % End if
        
        % Test Proximal linear
        [proxlinsol, proxlininfo] = proxlin(A, b, sqrt(maxiter * m), 0, init_x, ...
            maxiter, tol, true, 0, stepsize, 0, show_info);
        
        if proxlininfo.status == "Optimal"
            nProxLinEpochtoOpt(i, k) = proxlininfo.niter;
        end % End if 
        
        % Test Proximal point
        [proxptsol, proxptinfo] = proxpt(A, b, sqrt(maxiter * m), 0, init_x, ...
            maxiter, tol, true, 0, stepsize, 0, show_info);
        
        if proxptinfo.status == "Optimal"
            nProxPtEpochtoOpt(i, k) = proxptinfo.niter;
        end % End if 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Momentum
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test SGD
        [sgdsolm, sgdinfom] = proxsgd(A, b, sqrt(maxiter * m), beta, init_x, ...
            maxiter, tol, true, 1, 0, stepsize, 0, show_info);
        
        if sgdinfom.status == "Optimal"
            nSgdmEpochtoOpt(i, k) = sgdinfom.niter;
        end % End if
        
        % Test Proximal linear
        [proxlinsolm, proxlininfom] = proxlin(A, b, sqrt(maxiter * m), beta, init_x, ...
            maxiter, tol, true, 0, stepsize, 0, show_info);
        
        if proxlininfom.status == "Optimal"
            nProxLinmEpochtoOpt(i, k) = proxlininfom.niter;
        end % End if 
        
        % Test Proximal point
        [proxptsolm, proxptinfom] = proxpt(A, b, sqrt(maxiter * m), beta, init_x, ...
            maxiter, tol, true, 0, stepsize, 0, show_info);
        
        if proxptinfom.status == "Optimal"
            nProxPtmEpochtoOpt(i, k) = proxptinfom.niter;
        end % End if 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Minibatch
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test SGD
        [sgdsolb, sgdinfob] = proxsgd(A, b, sqrt(maxiter * m / batchsize), 0, init_x, ...
            maxiter, tol, true, batchsize, 0, stepsize, 0, show_info);
        
        if sgdinfob.status == "Optimal"
            nSgdbEpochtoOpt(i, k) = sgdinfob.niter;
        end % End if
        
        % Test Proximal linear
        [proxlinsolb, proxlininfob] = proxlinbatch(A, b, sqrt(maxiter * m / batchsize), 0, init_x, ...
            maxiter, tol, true, batchsize, 0, stepsize, show_info);
        
        if proxlininfob.status == "Optimal"
            nProxLinbEpochtoOpt(i, k) = proxlininfob.niter;
        end % End if 
        
        % Test Proximal point
        [proxptsolb, proxptinfob] = proxptbatch(A, b, sqrt(maxiter * m / batchsize), init_x, ...
            maxiter, tol, true, batchsize, 0, stepsize, show_info);
        
        if proxptinfob.status == "Optimal"
            nProxPtbEpochtoOpt(i, k) = proxptinfob.niter;
        end % End if 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Momentum and Minibatch
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test SGD
        [sgdsolbm, sgdinfobm] = proxsgd(A, b, sqrt(maxiter * m / batchsize), beta, init_x, ...
            maxiter, tol, true, batchsize, 0, stepsize, 0, show_info);
        
        if sgdinfobm.status == "Optimal"
            nSgdbmEpochtoOpt(i, k) = sgdinfobm.niter;
        end % End if
        
        % Test Proximal linear
        [proxlinsolbm, proxlininfobm] = proxlinbatch(A, b, sqrt(maxiter * m / batchsize), beta, init_x, ...
            maxiter, tol, true, batchsize, 0, stepsize, show_info);
        
        if proxlininfobm.status == "Optimal"
            nProxLinbmEpochtoOpt(i, k) = proxlininfobm.niter;
        end % End if 
        
        % Test Proximal point
        [proxptsolbm, proxptinfobm] = proxptbatch(A, b, sqrt(maxiter * m / batchsize), init_x, ...
            maxiter, tol, true, batchsize, 0, stepsize, show_info, beta);
        
        if proxptinfobm.status == "Optimal"
            nProxPtbmEpochtoOpt(i, k) = proxptinfobm.niter;
        end % End if 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end % End for
    
    fprintf("- Stepsize " + k + " done. \n"); 
    
end % End for

% Plot the summary graph
num_iter = maxiter * m;
xcord = steprange;

save("kappa_" + kappa + "_pfail_" + pfail + "_batch_" + batchsize + "_momentum_" + beta * 100 + "_epoch_env.mat");

hold off;
close all;

shadedErrorBarSemi(xcord, mean(nSgdEpochtoOpt / 300), std(nSgdEpochtoOpt / 300), "lineProps", {"-+", "MarkerSize", 18,  "LineWidth", 2});
% semilogx(xcord, sum(nSgdEpochtoOpt / (nTest * 300), 1), "-+", "MarkerSize", 18,  "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinEpochtoOpt / 300), std(nProxLinEpochtoOpt / 300), "lineProps", {"-o", "MarkerSize", 18,  "LineWidth", 2});
% semilogx(xcord, (sum(nProxLinEpochtoOpt / (nTest * 300), 1)), "-o", "MarkerSize", 18,  "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtEpochtoOpt / 300), std(nProxPtEpochtoOpt / 300), "lineProps", {"-d", "MarkerSize", 18,  "LineWidth", 2});
hold on;

shadedErrorBarSemi(xcord, mean(nSgdmEpochtoOpt / 300), std(nSgdmEpochtoOpt / 300), "lineProps", {"-s", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nSgdmEpochtoOpt / (nTest * 300), 1)), "-s", "MarkerSize", 18,  "LineWidth", 2 , "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinmEpochtoOpt / 300), std(nProxLinmEpochtoOpt / 300), "lineProps", {"-*", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nProxLinmEpochtoOpt / (nTest * 300), 1)), "-*", "MarkerSize", 18, "LineWidth", 2, "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtmEpochtoOpt / 300), std(nProxPtmEpochtoOpt / 300), "lineProps", {"-x", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});

set(gca, "FontSize", 20, "FontWeight", "bold")
xlim([min(steprange), max(steprange)]);

legend(["SGD", "SPL", "SPP", "SEGD", "SEPL", "SEPP"], "FontSize", 20);

savefig("kappa_" + kappa + "_pfail_" + pfail + "_batch_" + 1 + "_momentum_" + beta + "_epoch.fig");

hold off;
close all;

shadedErrorBarSemi(xcord, mean(nSgdbEpochtoOpt / 300), std(nSgdbEpochtoOpt / 300), "lineProps", {"-+", "MarkerSize", 18,  "LineWidth", 2});
% semilogx(xcord, sum(nSgdbEpochtoOpt / (nTest * 300), 1), "-+", "MarkerSize", 18,  "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinbEpochtoOpt / 300), std(nProxLinbEpochtoOpt / 300), "lineProps", {"-o", "MarkerSize", 18,  "LineWidth", 2});
% semilogx(xcord, (sum(nProxLinbEpochtoOpt / (nTest * 300), 1)), "-o", "MarkerSize", 18,  "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtbEpochtoOpt / 300), std(nProxPtbEpochtoOpt / 300), "lineProps", {"-d", "MarkerSize", 18,  "LineWidth", 2});
hold on;

shadedErrorBarSemi(xcord, mean(nSgdbmEpochtoOpt / 300), std(nSgdbmEpochtoOpt / 300), "lineProps", {"-s", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nSgdbmEpochtoOpt / (nTest * 300), 1)), "-s", "MarkerSize", 18,  "LineWidth", 2 , "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinbmEpochtoOpt / 300), std(nProxLinbmEpochtoOpt / 300), "lineProps", {"-*", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nProxLinbmEpochtoOpt / (nTest * 300), 1)), "-*", "MarkerSize", 18, "LineWidth", 2, "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtbmEpochtoOpt / 300), std(nProxPtbmEpochtoOpt / 300), "lineProps", {"-x", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});

set(gca, "FontSize", 20, "FontWeight", "bold")
xlim([min(steprange), max(steprange)]);

legend(["SGD", "SPL", "SPP", "SEGD", "SEPL", "SEPP"], "FontSize", 20);

savefig("kappa_" + kappa + "_pfail_" + pfail + "_batch_" + batchsize + "_momentum_" + beta + "_epoch.fig");


