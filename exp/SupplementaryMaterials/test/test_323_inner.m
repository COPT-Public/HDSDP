% Inner test script for Zipcode Phase Retrieval Momentum&Minibatch (3.2.3)
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
runb = true;

for k = 1:length(steprange)
    
    stepsize = steprange(k);
    runb = true;
    
    for i = 1:nTest
        
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
        if k >= 2
            [proxptsolb, proxptinfob] = proxptbatch(A, b, sqrt(maxiter * m / batchsize), init_x, ...
                maxiter, tol, true, batchsize, 0, stepsize, show_info);     
            if proxptinfob.status == "Optimal"
                nProxPtbEpochtoOpt(i, k) = proxptinfob.niter;
            end % End if
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

close all;

% Plot the summary graph
num_iter = maxiter * m;
xcord = steprange;

save("zipcode_" + idx + "_pfail_" + pfail + "_batch_" + batchsize + "_momentum_" + beta * 100 + "_epoch_env.mat");

hold off;
close all;

shadedErrorBarSemi(xcord, mean(nSgdEpochtoOpt / m), std(nSgdEpochtoOpt / m), "lineProps", {"-+", "MarkerSize", 18,  "LineWidth", 2});
% semilogx(xcord, sum(nSgdEpochtoOpt / (nTest * m), 1), "-+", "MarkerSize", 18,  "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinEpochtoOpt / m), std(nProxLinEpochtoOpt / m), "lineProps", {"-o", "MarkerSize", 18,  "LineWidth", 2});
% semilogx(xcord, (sum(nProxLinEpochtoOpt / (nTest * m), 1)), "-o", "MarkerSize", 18,  "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtEpochtoOpt / m), std(nProxPtEpochtoOpt / m), "lineProps", {"-d", "MarkerSize", 18,  "LineWidth", 2});
hold on;

shadedErrorBarSemi(xcord, mean(nSgdmEpochtoOpt / m), std(nSgdmEpochtoOpt / m), "lineProps", {"-s", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nSgdmEpochtoOpt / (nTest * m), 1)), "-s", "MarkerSize", 18,  "LineWidth", 2 , "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinmEpochtoOpt / m), std(nProxLinmEpochtoOpt / m), "lineProps", {"-*", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nProxLinmEpochtoOpt / (nTest * m), 1)), "-*", "MarkerSize", 18, "LineWidth", 2, "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtmEpochtoOpt / m), std(nProxPtmEpochtoOpt / m), "lineProps", {"-x", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});

set(gca, "FontSize", 20, "FontWeight", "bold")
xlim([min(steprange), max(steprange)]);
ylim([0, 400]);

legend(["SGD", "SPL", "SPP", "SEGD", "SEPL", "SEPP"], "FontSize", 20);
savefig("zipcode_" + idx + "_pfail_" + pfail + "_batch_" + 1 + "_momentum_" + beta + "_epoch.fig");

hold off;
close all;

shadedErrorBarSemi(xcord, mean(nSgdbEpochtoOpt / m), std(nSgdbEpochtoOpt / m), "lineProps", {"-+", "MarkerSize", 18,  "LineWidth", 2});
% semilogx(xcord, sum(nSgdbEpochtoOpt / (nTest * m), 1), "-+", "MarkerSize", 18,  "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinbEpochtoOpt / m), std(nProxLinbEpochtoOpt / m), "lineProps", {"-o", "MarkerSize", 18,  "LineWidth", 2});
% semilogx(xcord, (sum(nProxLinbEpochtoOpt / (nTest * m), 1)), "-o", "MarkerSize", 18,  "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtbEpochtoOpt / m), std(nProxPtbEpochtoOpt / m), "lineProps", {"-d", "MarkerSize", 18,  "LineWidth", 2});
hold on;

shadedErrorBarSemi(xcord, mean(nSgdbmEpochtoOpt / m), std(nSgdbmEpochtoOpt / m), "lineProps", {"-s", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nSgdbmEpochtoOpt / (nTest * m), 1)), "-s", "MarkerSize", 18,  "LineWidth", 2 , "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinbmEpochtoOpt / m), std(nProxLinbmEpochtoOpt / m), "lineProps", {"-*", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nProxLinbmEpochtoOpt / (nTest * m), 1)), "-*", "MarkerSize", 18, "LineWidth", 2, "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtbmEpochtoOpt / m), std(nProxPtbmEpochtoOpt / m), "lineProps", {"-x", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--"});

set(gca, "FontSize", 20, "FontWeight", "bold")
xlim([min(steprange), max(steprange)]);
ylim([0, 400]);

legend(["SGD", "SPL", "SPP", "SEGD", "SEPL", "SEPP"], "FontSize", 20);

savefig("zipcode_" + idx + "_pfail_" + pfail + "_batch_" + batchsize + "_momentum_" + beta + "_epoch.fig");


