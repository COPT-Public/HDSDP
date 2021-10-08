% Inner test script for Synthetic Blind Deconvolution Momentum&Minibatch (3.3.3)
U = data.U;
V = data.V;
b = data.b;
bestloss = data.bestloss;

tol = bestloss * 1.5;

[m, n] = size(U);

nSgdEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nSgdmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nSgdbEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nSgdbmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxLinEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxLinmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxLinbEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;
nProxLinbmEpochtoOpt = ones(nTest, length(steprange)) * maxiter * m;

for k = 1:length(steprange)
    
    stepsize = steprange(k);
    
    for i = 1:nTest
        
        init_z = randn(2 * n, 1);
        
        % Nothing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test SGD
        [sgdsol, sgdinfo] = proxsgdblind(U, V, b, sqrt(maxiter * m), 0, init_z, ...
            maxiter, tol, true, 0, stepsize, show_info);
        
        if sgdinfo.status == "Optimal"
            nSgdEpochtoOpt(i, k) = sgdinfo.niter;
        end % End if
        
        % Test Proximal linear
        [proxlinsol, proxlininfo] = proxlinblind(U, V, b, sqrt(maxiter * m), 0, init_z, maxiter, tol, ...
            early_stop, 0, stepsize, show_info);
        
        if proxlininfo.status == "Optimal"
            nProxLinEpochtoOpt(i, k) = proxlininfo.niter;
        end % End if 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Momentum
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test SGD
        [sgdsolm, sgdinfom] = proxsgdblind(U, V, b, sqrt(maxiter * m), beta, init_z, ...
            maxiter, tol, true, 0, stepsize, show_info);
        
        if sgdinfom.status == "Optimal"
            nSgdmEpochtoOpt(i, k) = sgdinfo.niter;
        end % End if
        
        % Test Proximal linear
        [proxlinsolm, proxlininfom] = proxlinblind(U, V, b, sqrt(maxiter * m), beta, init_z, maxiter, tol, ...
            early_stop, 0, stepsize, show_info);
        
        if proxlininfom.status == "Optimal"
            nProxLinmEpochtoOpt(i, k) = proxlininfo.niter;
        end % End if 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Minibatch
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test SGD
        [sgdsolb, sgdinfob] = proxsgdblindbatch(U, V, b, sqrt(maxiter * m / batchsize), 0, init_z, maxiter, tol, ...
                early_stop, batchsize, 0, stepsize, show_info);
    
        if sgdinfob.status == "Optimal"
            nSgdbEpochtoOpt(i, k) = sgdinfob.niter;
        end % End if
        
        % Test Proximal linear
        [proxlinsolb, proxlininfob] = proxlinblindbatch(U, V, b, sqrt(maxiter * m / batchsize), 0, init_z, maxiter, tol, ...
                early_stop, batchsize, 0, stepsize, show_info);
        
        if proxlininfob.status == "Optimal"
            nProxLinbEpochtoOpt(i, k) = proxlininfob.niter;
        end % End if 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Momentum and Minibatch
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test SGD
        [sgdsolbm, sgdinfobm] = proxsgdblindbatch(U, V, b, sqrt(maxiter * m / batchsize), beta, init_z, maxiter, tol, ...
                early_stop, batchsize, 0, stepsize, show_info);
        
        if sgdinfobm.status == "Optimal"
            nSgdbmEpochtoOpt(i, k) = sgdinfobm.niter;
        end % End if
        
        % Test Proximal linear
        [proxlinsolbm, proxlininfobm] = proxlinblindbatch(U, V, b, sqrt(maxiter * m / batchsize), beta, init_z, maxiter, tol, ...
                early_stop, batchsize, 0, stepsize, show_info);
        
        if proxlininfobm.status == "Optimal"
            nProxLinbmEpochtoOpt(i, k) = proxlininfobm.niter;
        end % End if 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end % End for
    
    fprintf("- Stepsize " + k + " done. \n"); 
    
end % End for

% Plot the summary graph
num_iter = maxiter * m;
xcord = steprange;

semilogx(xcord, sum(nSgdEpochtoOpt / (nTest * 300), 1), "-+", "MarkerSize", 18,  "LineWidth", 2);
hold on;

semilogx(xcord, (sum(nProxLinEpochtoOpt / (nTest * 300), 1)), "-o", "MarkerSize", 18,  "LineWidth", 2);
hold on;

semilogx(xcord, (sum(nSgdmEpochtoOpt / (nTest * 300), 1)), "-s", "MarkerSize", 18,  "LineWidth", 2 , "LineStyle", "--");
hold on;

semilogx(xcord, (sum(nProxLinmEpochtoOpt / (nTest * 300), 1)), "-*", "MarkerSize", 18,  "LineWidth", 2, "LineStyle", "--");
hold on;

set(gca, "FontSize", 20, "FontWeight", "bold")
xlim([min(steprange), max(steprange)]);

legend(["SGD", "SPL", "SEGD", "SEPL"], "FontSize", 20);

savefig("blind_kappa_" + kappa + "_pfail_" + pfail + "_batch_" + 1 + "_momentum_" + beta + "_epoch_env.fig");

hold off;
close all;

semilogx(xcord, sum(nSgdbEpochtoOpt / (nTest * 300), 1), "-+", "MarkerSize", 18,  "LineWidth", 2);
hold on;

semilogx(xcord, (sum(nProxLinbEpochtoOpt / (nTest * 300), 1)), "-o", "MarkerSize", 18,  "LineWidth", 2);
hold on;

semilogx(xcord, (sum(nSgdbmEpochtoOpt / (nTest * 300), 1)), "-s", "MarkerSize", 18,  "LineWidth", 2 , "LineStyle", "--");
hold on;

semilogx(xcord, (sum(nProxLinbmEpochtoOpt / (nTest * 300), 1)), "-*", "MarkerSize", 18, "LineWidth", 2, "LineStyle", "--");
hold on;

set(gca, "FontSize", 20, "FontWeight", "bold")
xlim([min(steprange), max(steprange)]);

legend(["SGD", "SPL", "SEGD", "SEPL"], "FontSize", 20);

savefig("blind_kappa_" + kappa + "_pfail_" + pfail + "_batch_" + batchsize + "_momentum_" + beta + "_epoch.fig");
hold off;
