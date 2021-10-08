% Outer test script for Synthetic Blind Deconvolution Minibatch&Momentum (3.3.3)
clear;
clc;
close all;

addpath(fullfile("..", "opt"));
addpath(fullfile("..", "data"));

% Set range of batchsizes of interest (excluding 1)
batchrange = [32];

% Set range of initial stepsizes of interest
% e.g., from 0.1 to 100 under logscale
steprange = logspace(-2, 1, 10);
nsteptotest = length(steprange);

% Momentum Range
momentumrange = [0.2, 0.6];

% Set range of kappa values of interest
kapparange = [10];

% Set range of pfail values of interest (here we use pfail * 10 for
% representation)
pfailrange = [2, 3];

% Set number of repeats
nTest = 20;

% Set optimization parameters
use_vm = 0;
maxiter = 400;
show_info = false;
early_stop = true;

% Start experiment
% Set global random seed
rng(123);

for pfail = pfailrange
    for kappa = kapparange
        for batchsize = batchrange
            for beta = momentumrange
                close all;
                
                fprintf("Start experiment for pfail = " + pfail / 10 + ", kappa = " + ...
                    kappa + ", batchsize = " + batchsize + ", momentum = " + beta + "\n");
                
                file = fullfile("..", "data", "blind_kappa_" + kappa + "_pfail_" + pfail + ".mat");
                load(file);
                
                test_333_inner;
                
                fprintf("End experiment for pfail = " + pfail / 10 + ", kappa = " + ...
                    kappa + ", batchsize = " + batchsize + ", momentum = " + beta + "\n");
                
            end % End for
        end % End for
    end % End for
end % End for


