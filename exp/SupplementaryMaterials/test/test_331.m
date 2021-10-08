% Test script for Synthetic Blind Deconvolution Momentum (3.3.1)
clear;
clc;
close all;

addpath(fullfile("..", "opt"));
addpath(fullfile("..", "data"));

kappa = 1;
pfail = 0.3;

file = "blind_kappa_" + kappa + "_pfail_" + pfail * 10 + ".mat";
load(file);

alpharange = logspace(-2, 0, 10);
nTest = 20;

U = data.U;
V = data.V;
b = data.b;

maxiter = 400;
obj_beta = 0.2 ;
init_z = randn(200, 1);
gamma = sqrt(maxiter * 300);
tol = data.bestloss * 1.5;
early_stop = true;
adp_param = 0;
alpha_0 = 1;

nSgdEpochtoOpt = zeros(nTest, length(alpharange));
nSgdmEpochtoOpt = zeros(nTest, length(alpharange));
nProxLinEpochtoOpt = zeros(nTest, length(alpharange));
nProxLinmEpochtoOpt = zeros(nTest, length(alpharange));
nProxPtEpochtoOpt = zeros(nTest, length(alpharange));
nProxPtmEpochtoOpt = zeros(nTest, length(alpharange));

show_info = false;

for k = 1:length(alpharange)
    
    alpha_0 = alpharange(k);
    
    parfor i = 1:nTest
        
        beta = 0;
        
        [sgdsol, sgdinfo] = proxsgdblind(U, V, b, gamma, beta, init_z, maxiter, tol, ...
            early_stop, adp_param, alpha_0, show_info);
        
        nSgdEpochtoOpt(i, k) = sgdinfo.nepochs;
        
        [proxlinsol, proxlininfo] = proxlinblind(U, V, b, gamma, beta, init_z, maxiter, tol, ...
            early_stop, adp_param, alpha_0, show_info);
        
        nProxLinEpochtoOpt(i, k) = proxlininfo.nepochs;
        
        [proxptsol, proxptinfo] = proxptblind(U, V, b, gamma, beta, init_z, maxiter, tol, ...
            early_stop, adp_param, alpha_0, show_info);
        
        nProxPtEpochtoOpt(i, k) = proxptinfo.nepochs;
        
        beta = obj_beta;
        
        [sgdsolm, sgdinfom] = proxsgdblind(U, V, b, gamma, beta, init_z, maxiter, tol, ...
            early_stop, adp_param, alpha_0, show_info);
        
        nSgdmEpochtoOpt(i, k) = sgdinfom.nepochs;
        
        [proxlinsolm, proxlininfom] = proxlinblind(U, V, b, gamma, beta, init_z, maxiter, tol, ...
            early_stop, adp_param, alpha_0, show_info);
        
        nProxLinmEpochtoOpt(i, k) = proxlininfom.nepochs;
        
        [proxptsolm, proxptinfom] = proxptblind(U, V, b, gamma, beta, init_z, maxiter, tol, ...
            early_stop, adp_param, alpha_0, show_info);
        
        nProxPtmEpochtoOpt(i, k) = proxptinfom.nepochs;
        
    end % End for
    
    disp("Done once");
    
end % End for

% Plot the summary graph
[m, ~] = size(data.U);
num_iter = maxiter * m;
xcord = alpharange;

semilogx(xcord, sum(nSgdEpochtoOpt / nTest, 1), "-+", "LineWidth", 2);
hold on;

semilogx(xcord, (sum(nProxLinEpochtoOpt / nTest, 1)), "-o", "LineWidth", 2);
hold on;

semilogx(xcord, (sum(nProxPtEpochtoOpt / nTest, 1)), "-x", "LineWidth", 2);
hold on;

semilogx(xcord, (sum(nSgdmEpochtoOpt / nTest, 1)), "-s", "LineWidth", 2 , "LineStyle", "--");
hold on;

semilogx(xcord, (sum(nProxLinmEpochtoOpt / nTest, 1)), "-*", "LineWidth", 2, "LineStyle", "--");
hold on;

semilogx(xcord, (sum(nProxPtmEpochtoOpt / nTest, 1)), "-d", "LineWidth", 2, "LineStyle", "--");
hold on;

set(gca, "FontSize", 20, "FontWeight", "bold")
xlim([min(alpharange), max(alpharange)]);

% xlabel("alpha0");
% ylabel("epochs");

legend(["SGD", "SPL", "SPP", "SEGD",...
    "SEPL", "SEPP"], "FontSize", 20);

save("blind_kappa_" + kappa + "_pfail_" + pfail + "_momentum_" + obj_beta + "_epoch.mat");
savefig("blind_kappa_" + kappa + "_pfail_" + pfail + "_momentum_" + obj_beta * 100 + "_epoch.fig");

hold off;
