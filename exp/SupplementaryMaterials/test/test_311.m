% Test script for Synthetic Phase Retrieval Momentum (3.1.1)
% Code for camera-ready paper: Passed by Gwz
% This experiment corresponds to the first two sub-figures
% in Figure 3 of the paper

close all;
clear;
clc;

addpath(fullfile("..", "opt"));
addpath(fullfile("..", "data"));

kappa = 10;
pfail = 2;
beta = 0.6;
file = fullfile("..", "data", "kappa_" + kappa + "_pfail_" + pfail + ".mat");
nTest = 20;

alpharange = logspace(-2, 0, 10);

load(file)

bestloss = data.bestloss;
params = setdefaultparams(data);
params.use_vm = 0;
params.beta = 0;
params.maxiter = 400;
params.gamma = 1/sqrt(params.maxiter * 300);
params.adp_param = 0;
params.alpha_0 = 1;
params.early_stop = true;
params.logscale = false;
params.show_info = false;

nSgdEpochtoOpt = zeros(nTest, length(alpharange));
nSgdmEpochtoOpt = zeros(nTest, length(alpharange));
nProxLinEpochtoOpt = zeros(nTest, length(alpharange));
nProxLinmEpochtoOpt = zeros(nTest, length(alpharange));
nProxPtEpochtoOpt = zeros(nTest, length(alpharange));
nProxPtmEpochtoOpt = zeros(nTest, length(alpharange));

for k = 1:length(alpharange)
    params.alpha_0 = alpharange(k);
    for i = 1:nTest
        params.init_x = randn(100, 1);
        params.beta = 0;
        % Test SGD
        params.optmethod = "SGD";
        [sgdsol, sgdinfo] = proxopt(data, params);
        nSgdEpochtoOpt(i, k) = sgdinfo.nepochs;
        
        % Test Proximal linear
        params.optmethod = "ProxLin";
        [proxlinsol, proxlininfo] = proxopt(data, params);
        
        nProxLinEpochtoOpt(i, k) = proxlininfo.nepochs;
        
        % Test Proximal point
        params.optmethod = "ProxPt";
        [proxptsol, proxptinfo] = proxopt(data, params);
        nProxPtEpochtoOpt(i, k) = proxptinfo.nepochs;
        
        % Test methods with momentum
        params.beta = beta;
        
        % Test SGD
        params.optmethod = "SGD";
        [sgdsolm, sgdinfom] = proxopt(data, params);
        nSgdmEpochtoOpt(i, k) = sgdinfom.nepochs;
        
        % Test Proximal linear
        params.optmethod = "ProxLin";
        [proxlinsolm, proxlininfom] = proxopt(data, params);
        nProxLinmEpochtoOpt(i, k) = proxlininfom.nepochs;
        
        % Test Proximal point
        params.optmethod = "ProxPt";
        [proxptsolm, proxptinfom] = proxopt(data, params);
        nProxPtmEpochtoOpt(i, k) = proxptinfom.nepochs;
        
    end % End for
end % End for

% Plot the summary graph
[m, ~] = size(data.A);
xcord = alpharange;

shadedErrorBarSemi(xcord, mean(nSgdEpochtoOpt), std(nSgdEpochtoOpt), "lineProps", {"-+", "LineWidth", 2});
% semilogx(xcord, sum(nSgdEpochtoOpt / nTest, 1), "-+", "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinEpochtoOpt), std(nProxLinEpochtoOpt), "lineProps", {"-o", "LineWidth", 2});
% semilogx(xcord, (sum(nProxLinEpochtoOpt / nTest, 1)), "-o", "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtEpochtoOpt), std(nProxPtEpochtoOpt), "lineProps", {"-x", "LineWidth", 2});
% semilogx(xcord, (sum(nProxPtEpochtoOpt / nTest, 1)), "-x", "LineWidth", 2);
hold on;

shadedErrorBarSemi(xcord, mean(nSgdmEpochtoOpt), std(nSgdmEpochtoOpt), "lineProps", {"-s", "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nSgdmEpochtoOpt / nTest, 1)), "-s", "LineWidth", 2 , "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxLinmEpochtoOpt), std(nProxLinmEpochtoOpt), "lineProps", {"-*", "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nProxLinmEpochtoOpt / nTest, 1)), "-*", "LineWidth", 2, "LineStyle", "--");
hold on;

shadedErrorBarSemi(xcord, mean(nProxPtmEpochtoOpt), std(nProxPtmEpochtoOpt), "lineProps", {"-d", "LineWidth", 2, "LineStyle", "--"});
% semilogx(xcord, (sum(nProxPtmEpochtoOpt / nTest, 1)), "-d", "LineWidth", 2, "LineStyle", "--");
hold on;
set(gca, "FontSize", 20, "FontWeight", "bold")
xlim([min(alpharange), max(alpharange)]);

% xlabel("alpha0");
% ylabel("epochs");

legend(["SGD", "SPL", "SPP", "SEGD", "SEPL", "SPP"], "FontSize", 20);
tightfig;
ylim([0, 400]);

savefig("kappa_" + kappa + "_pfail_" + pfail + "_momentum_" + beta * 100 + "_epoch");
close all;
save("kappa_" + kappa + "_pfail_" + pfail + "_momentum_" + beta * 100 + "_epoch.mat");

hold off;
