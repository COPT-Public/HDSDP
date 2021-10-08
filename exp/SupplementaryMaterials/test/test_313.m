% Outer test script for Synthetic Phase Retrieval Momentum&Minibatch (3.1.3)
% The first two sub-figures from Figure 4
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
momentumrange = [0.6];

% Set range of kappa values of interest
kapparange = [10];

% Set range of pfail values of interest (here we use pfail * 10 for
% representation)
pfailrange = [3];

% Set number of repeats
nTest = 20;

% Set optimization parameters
use_vm = 0;
maxiter = 400;
show_info = false;

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
                
                file = fullfile("..", "data", "kappa_" + kappa + "_pfail_" + pfail + ".mat");
                load(file);
                
                test_313_inner;
                
                fprintf("End experiment for pfail = " + pfail / 10 + ", kappa = " + ...
                    kappa + ", batchsize = " + batchsize + ", momentum = " + beta + "\n");
                
            end % End for
        end % End for
    end % End for
end % End for


