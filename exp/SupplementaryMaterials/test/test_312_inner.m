% Inner test script for Synthetic Phase Retrieval Minibatch (3.1.2)
% Experiment setup should be isolated in console file
% Variables already provided in outer environment
% 1. pfail 2. kappa 3. batchrange 4. steprange 5. data
A = data.A;
b = data.b;
bestloss = data.bestloss;

tol = bestloss * 1.5;

[m, n] = size(A);

% Initialize result arrays with maximum number of iterations
nsgdIter = ones(nbatchtotest + 1, nsteptotest, nTest);
nproxlinIter = ones(nbatchtotest + 1, nsteptotest, nTest);
nproxptIter = ones(nbatchtotest + 1, nsteptotest, nTest);
nsgdIterToOpt = ones(nbatchtotest + 1, nsteptotest) * maxiter * m;
nproxlinIterToOpt = ones(nbatchtotest + 1, nsteptotest) * maxiter * m;
nproxptIterToOpt = ones(nbatchtotest + 1, nsteptotest) * maxiter * m;
nsgdIterToOptStd = ones(nbatchtotest + 1, nsteptotest) * maxiter * m;
nproxlinIterToOptStd = ones(nbatchtotest + 1, nsteptotest) * maxiter * m;
nproxptIterToOptStd = ones(nbatchtotest + 1, nsteptotest) * maxiter * m;

% Do experiment with batchsize 1
batchsize = 1;

idx = 0;
for stepsize = steprange
    
    idx = idx + 1;
    tempSgdIter = ones(nTest, 1) * maxiter * m;
    tempProxlinIter = ones(nTest, 1) * maxiter * m;
    tempProxptIter = ones(nTest, 1) * maxiter * m;
    
    parfor i = 1:nTest
        
        init_x = randn(n, 1);
        
        [sgdsol, sgdinfo] = proxsgd(A, b, sqrt(maxiter * m), 0, init_x, ...
            maxiter, tol, true, batchsize, 0, stepsize, 0, show_info);
        
        [proxlinsol, proxlininfo] = proxlin(A, b, sqrt(maxiter * m), 0, init_x, ...
            maxiter, tol, true, 0, stepsize, 0, show_info);
        
        [proxptsol, proxptinfo] = proxpt(A, b, sqrt(maxiter * m), 0, init_x, ...
            maxiter, tol, true, 0, stepsize, 0, show_info);
        
        if sgdinfo.status == "Optimal"
            tempSgdIter(i) = sgdinfo.niter;
        end % End if
        
        if proxlininfo.status == "Optimal"
            tempProxlinIter(i) = proxlininfo.niter;
        end % End if
        
        if proxptinfo.status == "Optimal"
            tempProxptIter(i) = proxptinfo.niter;
        end % End if
        
    end % End parfor
    
    nsgdIterToOpt(1, idx) = mean(tempSgdIter);
    nproxlinIterToOpt(1, idx) = mean(tempProxlinIter);
    nproxptIterToOpt(1, idx) = mean(tempProxptIter);
    nsgdIterToOptStd(1, idx) = std(tempSgdIter);
    nproxlinIterToOptStd(1, idx) = std(tempProxlinIter);
    nproxptIterToOptStd(1, idx) = std(tempProxptIter);
    nsgdIter(1, idx, :) = tempSgdIter;
    nproxlinIter(1, idx, :) = tempProxlinIter;
    nproxptIter(1, idx, :) = tempProxptIter;
    
end % End for

fprintf("- Batchsize 1 done \n");

for k = 1:nbatchtotest
    
    batchsize = batchrange(k);
    
    idx = 0;
    for stepsize = steprange
        
        idx = idx + 1;
        tempSgdIter = ones(nTest, 1) * maxiter * m;
        tempProxlinIter = ones(nTest, 1) * maxiter * m;
        tempProxptIter = ones(nTest, 1) * maxiter * m;
        
        parfor i = 1:nTest
            
            init_x = randn(n, 1);
            
            [sgdsol, sgdinfo] = proxsgd(A, b, sqrt(maxiter * m / batchsize), 0, init_x, ...
                maxiter, tol, true, batchsize, 0, stepsize, 0, show_info);
            
            [proxlinsol, proxlininfo] = proxlinbatch(A, b, sqrt(maxiter * m / batchsize), 0, init_x, ...
                maxiter, tol, true, batchsize, 0, stepsize, show_info);
            
            [proxptsol, proxptinfo] = proxptbatch(A, b, sqrt(maxiter * m / batchsize), init_x, ...
                maxiter, tol, true, batchsize, 0, stepsize, show_info);
            
            if sgdinfo.status == "Optimal"
                tempSgdIter(i) = sgdinfo.niter;
            end % End if
            
            if proxlininfo.status == "Optimal"
                tempProxlinIter(i) = proxlininfo.niter;
            end % End if
            
            if proxptinfo.status == "Optimal"
                tempProxptIter(i) = proxptinfo.niter;
            end % End if
            
        end % End parfor
        
        nsgdIterToOpt(k + 1, idx) = mean(tempSgdIter);
        nproxlinIterToOpt(k + 1, idx) = mean(tempProxlinIter);
        nproxptIterToOpt(k + 1, idx) = mean(tempProxptIter);
        nsgdIterToOptStd(k + 1, idx) = std(tempSgdIter);
        nproxlinIterToOptStd(k + 1, idx) = std(tempProxlinIter);
        nproxptIterToOptStd(k + 1, idx) = std(tempProxptIter);
        nsgdIter(k + 1, idx, :) = tempSgdIter;
        nproxlinIter(k + 1, idx, :) = tempProxlinIter;
        nproxptIter(k + 1, idx, :) = tempProxptIter;
        
    end % End for
    
    fprintf("- Batchsize " + batchrange(k) + " done \n");
    
end % End for

envname = "kappa_" + kappa + "_pfail_" + pfail + "_env.mat";
fprintf("Saving environment to " + envname + ". \n");
save(envname);

sgdSpeedup = ones(nTest, nbatchtotest + 1);
nsgdBatchOneIter = ones(nTest, 1);
proxlinSpeedup = ones(nTest, nbatchtotest + 1);
nproxlinBatchOneIter = ones(nTest, 1);
proxptSpeedup = ones(nTest, nbatchtotest + 1);
nproxptBatchOneIter = ones(nTest, 1);

for i = 1:nTest
    tempSpeedup = min((nsgdIter(:, :, i)'));
    nsgdBatchOneIter(i) = tempSpeedup(1);
    sgdSpeedup(i, :) = max(tempSpeedup) ./ tempSpeedup;
    tempSpeedup = min((nproxlinIter(:, :, i)'));
    nproxlinBatchOneIter(i, :) = tempSpeedup(1);
    proxlinSpeedup(i, :) = max(tempSpeedup) ./ tempSpeedup;
    tempSpeedup = min((nproxptIter(:, :, i)'));
    nproxptBatchOneIter(i) = tempSpeedup(1);
    proxptSpeedup(i, :) = max(tempSpeedup) ./ tempSpeedup;
end % End for

fprintf("Experiments ended. Start plotting. \n");
fprintf("Type 1 figure: Speedup vs. Batchsize \n");

% Plot figure type 1: speedup
shadedErrorBar([1, batchrange], mean(sgdSpeedup), std(sgdSpeedup), "lineProps", {"-x",  "LineWidth", 3});
% plot([1, batchrange], proxlinSpeedup, "-x",  "LineWidth", 3);
hold on;
shadedErrorBar([1, batchrange], mean(proxlinSpeedup), std(proxlinSpeedup), "lineProps", {"-o",  "LineWidth", 3});
% plot([1, batchrange], sgdSpeedup, "-o", "LineWidth", 3);
hold on;
shadedErrorBar([1, batchrange], mean(proxptSpeedup), std(proxptSpeedup), "lineProps", {"-+",  "LineWidth", 3});
% plot([1, batchrange], proxptSpeedup, "-+", "LineWidth", 3);
legend(["SPL", "SGD", "SPP"], "FontSize", 20);
% xlabel("batchsize m");
% ylabel("speedup");
set(gca, "Fontsize", 20);
savefig("kappa_" + kappa + "_pfail_" + pfail + "_type_1.fig")
fprintf("Done. \n");

close all;

fprintf("Type 2 figure: Speedup vs. Stepsize \n");

% Plot figure type 2: robustness
for i = 1:nrobusttest
    
    idx = robustnessbatchidx(i) + 1;
    batchsize = robustnessbatch(i);
    
    fprintf("- Type 2 figure : Batchsize " + batchsize + " \n");
    robustsgdSpeedup = zeros(nTest, nsteptotest);
    robustproxlinSpeedup = zeros(nTest, nsteptotest);
    robustproxptSpeedup = zeros(nTest, nsteptotest);
    
    for j = 1:nTest
        robustsgdSpeedup(j, :) = nsgdBatchOneIter(j) ./ nsgdIter(idx, :, j);
        robustproxlinSpeedup(j, :) = nproxlinBatchOneIter(j) ./ nproxlinIter(idx, :, j);
        robustproxptSpeedup(j, :) = nproxptBatchOneIter(j) ./ nproxptIter(idx, :, j);
    end % End for
    
    % semilogx(steprange, robustproxlinSpeedup, "-x", "LineWidth", 3);
    shadedErrorBarSemi(steprange, mean(robustproxlinSpeedup), std(robustproxlinSpeedup), "lineProps", {"-x", "LineWidth", 3});
    hold on;
    % semilogx(steprange, robustsgdSpeedup, "-o", "LineWidth", 3);
    shadedErrorBarSemi(steprange, mean(robustsgdSpeedup), std(robustsgdSpeedup), "lineProps", {"-o", "LineWidth", 3});
    hold on;
    % semilogx(steprange, robustproxptSpeedup, "-+", "LineWidth", 3);
    shadedErrorBarSemi(steprange, mean(robustproxptSpeedup), std(robustproxptSpeedup), "lineProps", {"-+", "LineWidth", 3});
    
    legend(["SPL", "SGD", "SPP"], "FontSize", 20);
    set(gca, "Fontsize", 20);
    savefig("kappa_" + kappa + "_pfail_" + pfail + "_type_2_batch_" + batchsize + ".fig")
    
    %     xlabel("stepsize");
    %     ylabel("speedup");
    
    hold off;
    close all;
    
end % End

fprintf("Done. \n");

fprintf("Type 3 figure: Iteration to Opt vs. Stepsize");

% Plot figure type 3: Iteration to Opt
for i = 1:nrobusttest
    
    idx = robustnessbatchidx(i) + 1;
    batchsize = robustnessbatch(i);
    
    fprintf("- Type 3 figure : Batchsize " + batchsize + " \n");
    
    shadedErrorBarLoglog(steprange, nproxlinIterToOpt(idx, :), ...
        std(transpose(reshape(nproxlinIter(idx, :, :), nsteptotest, nTest))), "lineProps", {"-x", "LineWidth", 3});
    % semilogx(steprange, nproxlinIterToOpt(idx, :), "-x", "LineWidth", 3);
    hold on;
    shadedErrorBarLoglog(steprange, nsgdIterToOpt(idx, :), ...
        std(transpose(reshape(nsgdIter(idx, :, :), nsteptotest, nTest))), "lineProps", {"-x", "LineWidth", 3});
    % semilogx(steprange, nsgdIterToOpt(idx, :), "-o", "LineWidth", 3);
    hold on;
    shadedErrorBarLoglog(steprange, nproxptIterToOpt(idx, :), ...
        std(transpose(reshape(nproxptIter(idx, :, :), nsteptotest, nTest))), "lineProps", {"-x", "LineWidth", 3});
    % semilogx(steprange, nproxptIterToOpt(idx, :), "-+", "LineWidth", 3);
    
    legend(["SPL", "SGD", "SPP"], "FontSize", 20);
    set(gca, "Fontsize", 20);
    savefig("kappa_" + kappa + "_pfail_" + pfail + "_type_3_batch_" + batchsize + ".fig")
    
    %     xlabel("stepsize");
    %     ylabel("speedup");
    
    hold off;
    close all;
    
end % End

fprintf("Done. \n");

fprintf("**************************\n");
fprintf("*********  Done  *********\n");
fprintf("**************************\n");

