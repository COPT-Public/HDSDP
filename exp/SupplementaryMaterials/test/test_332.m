% Outer test script for Synthetic Blind Deconvolution Minibatch (3.3.2)
clear;
clc;
close all;

addpath(fullfile("..", "opt"));
addpath(fullfile("..", "data"));

% Set range of batchsizes of interest (excluding 1)
batchrange = [4, 8, 16, 32, 64];
robustnessbatch = [8, 32];
robustnessbatchidx = [2, 4];
nbatchtotest = length(batchrange);
nrobusttest = length(robustnessbatch);

% Set range of initial stepsizes of interest
% e.g., from 0.1 to 100 under logscale
steprange = logspace(-1, 2, 10);
nsteptotest = length(steprange);

% Set range of kappa values of interest
kapparange = [10];

% Set range of pfail values of interest (here we use pfail * 10 for 
% representation)
pfailrange = [2, 3];

% Set number of repeats
nTest = 20;

% Set optimization parameters
use_vm = 0;
beta = 0;
maxiter = 200;
show_info = false;
early_stop = true;

% Start experiment
% Set global random seed
rng(123);

for pfail = pfailrange
    for kappa = kapparange
        
        close all;

        fprintf("Start experiment for pfail = " + pfail / 10 + ", kappa = " + ...
            kappa + "\n");
        
        file = "blind_kappa_" + kappa + "_pfail_" + pfail + ".mat";
        load(file);

        test_332_inner;
        
        fprintf("Ended experiment for pfail = " + pfail / 10 + ", kappa = " + ...
            kappa + "\n");
    end % End for
end % End for


