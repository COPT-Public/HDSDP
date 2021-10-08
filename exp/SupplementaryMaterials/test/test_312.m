% Outer test script for Synthetic Phase Retrieval Minibatch (3.1.2)
% Code for camera-ready paper: Passed by Gwz
% This experiment corresponds to the first two and first four sub-figures
% in Figure 1 and Figure 2 of the paper

clear;
clc;
close all;

addpath(fullfile("..", "opt"));
addpath(fullfile("..", "data"));

% Set range of batchsizes of interest (excluding 1)
batchrange = [32, 64]; % [4, 8, 16, 32, 64];
robustnessbatch = [32, 64]; % [8, 32];
robustnessbatchidx = [1, 2]; % [2, 4];
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
nTest = 2;

% Set optimization parameters
use_vm = 0;
beta = 0;

% maxiter here denotes number of epochs and the number of inner iterations
% is given by m * maxiter
maxiter = 200;
show_info = false;

% Start experiment
% Set global random seed
rng(123);

for pfail = pfailrange
    for kappa = kapparange
        
        close all;

        fprintf("Start experiment for pfail = " + pfail / 10 + ", kappa = " + ...
            kappa + "\n");
        
        file = fullfile("..", "data", "kappa_" + kappa + "_pfail_" + pfail + ".mat");
        load(file);

        test_312_inner;
        
        fprintf("Ended experiment for pfail = " + pfail / 10 + ", kappa = " + ...
            kappa + "\n");
    end % End for
end % End for


