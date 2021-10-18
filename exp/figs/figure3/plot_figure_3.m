hold off;
close all;

shadedErrorBarSemi(xcord, mean(nSgdEpochtoOpt / 300), std(nSgdEpochtoOpt / 300), "lineProps", {"-+", "MarkerSize", 18,  "LineWidth", 2, "Color", [0.12, 0.56, 1.00]});
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

shadedErrorBarSemi(xcord, mean(nSgdbEpochtoOpt / 300), std(nSgdbEpochtoOpt / 300), "lineProps", {"-+", "MarkerSize", 18,  "LineWidth", 2, "Color", [0.12, 0.56, 1.00]});
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
