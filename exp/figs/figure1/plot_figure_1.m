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

fprintf("Type 1 figure: Speedup vs. Batchsize \n");

% Plot figure type 1: speedup
shadedErrorBar([1, batchrange], mean(sgdSpeedup), std(sgdSpeedup), "lineProps", {"-d",  "LineWidth", 3, "Color", [0.12, 0.56, 1.00]});
% plot([1, batchrange], proxlinSpeedup, "-x",  "LineWidth", 3);
hold on;
shadedErrorBar([1, batchrange], mean(proxlinSpeedup), std(proxlinSpeedup), "lineProps", {"-o",  "LineWidth", 3, "Color", [0.13, 0.55, 0.13]});
% plot([1, batchrange], sgdSpeedup, "-o", "LineWidth", 3);
hold on;
shadedErrorBar([1, batchrange], mean(proxptSpeedup), std(proxptSpeedup), "lineProps", {"-+",  "LineWidth", 3, "Color", "r"});
% plot([1, batchrange], proxptSpeedup, "-+", "LineWidth", 3);
legend(["SPL", "SGD", "SPP"], "FontSize", 20);
% xlabel("batchsize m");
% ylabel("speedup");
set(gca, "Fontsize", 20);
savefig("zipcode_" + imgidx + "_pfail_" + pfail + "_type_1.fig")
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
    
    for j = 1:16
        robustsgdSpeedup(j, :) = nsgdBatchOneIter(j) ./ nsgdIter(idx, :, j);
        robustproxlinSpeedup(j, :) = nproxlinBatchOneIter(j) ./ nproxlinIter(idx, :, j);
        if j <= 16
        robustproxptSpeedup(j, :) = nproxptBatchOneIter(j) ./ nproxptIter(idx, :, j);
        end % End if 
    end % End for
    
    % semilogx(steprange, robustproxlinSpeedup, "-x", "LineWidth", 3);
    shadedErrorBarSemi(steprange, mean(robustproxlinSpeedup), std(robustproxlinSpeedup), "lineProps", {"-d", "LineWidth", 3, "Color", [0.12, 0.56, 1.00]});
    hold on;
    % semilogx(steprange, robustsgdSpeedup, "-o", "LineWidth", 3);
    shadedErrorBarSemi(steprange, mean(robustsgdSpeedup), std(robustsgdSpeedup), "lineProps", {"-o", "LineWidth", 3, "Color", [0.13, 0.55, 0.13]});
    hold on;
    % semilogx(steprange, robustproxptSpeedup, "-+", "LineWidth", 3);
    shadedErrorBarSemi(steprange, mean(robustproxptSpeedup), std(robustproxptSpeedup), "lineProps", {"-+", "LineWidth", 3, "Color", "r"});
    
    legend(["SPL", "SGD", "SPP"], "FontSize", 20);
    set(gca, "Fontsize", 20);
    savefig("zipcode_" + imgidx + "_pfail_" + pfail + "_type_2_batch_" + batchsize + ".fig")
    
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
        std(transpose(reshape(nproxlinIter(idx, :, :), nsteptotest, nTest))), "lineProps", {"-d", "LineWidth", 3, "Color", [0.12, 0.56, 1.00]});
    % semilogx(steprange, nproxlinIterToOpt(idx, :), "-x", "LineWidth", 3);
    hold on;
    shadedErrorBarLoglog(steprange, nsgdIterToOpt(idx, :), ...
        min(std(transpose(reshape(nsgdIter(idx, :, :), nsteptotest, nTest))), nsgdIterToOpt(idx, :)), "lineProps", {"-x", "LineWidth", 3, "Color", [0.13, 0.55, 0.13]});
    % semilogx(steprange, nsgdIterToOpt(idx, :), "-o", "LineWidth", 3);
    hold on;
    shadedErrorBarLoglog(steprange, nproxptIterToOpt(idx, :), ...
        std(transpose(reshape(nproxptIter(idx, :, :), nsteptotest, nTest))), "lineProps", {"-x", "LineWidth", 3, "Color", "r"});
    % semilogx(steprange, nproxptIterToOpt(idx, :), "-+", "LineWidth", 3);
    
    legend(["SPL", "SGD", "SPP"], "FontSize", 20);
    set(gca, "Fontsize", 20);
    savefig("zipcode_" + imgidx + "_pfail_" + pfail + "_type_3_batch_" + batchsize + ".fig")
    
    %     xlabel("stepsize");
    %     ylabel("speedup");
    
    hold off;
    close all;
    
end % End

fprintf("Done. \n");