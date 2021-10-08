% Outer test script for Zipcode Phase Retrieval Momentum&Minibatch (3.2.3)
% The last two sub-figures from Figure 4
clear;
clc;
close all;

addpath(fullfile("..", "opt"));
addpath(fullfile("..", "data"));

% Set range of batchsizes of interest (excluding 1)
batchrange = [32];

% Set range of initial stepsizes of interest
% e.g., from 0.1 to 100 under logscale
steprange = logspace(0, 1, 10);
nsteptotest = length(steprange);

% Momentum Range
momentumrange = [0.9];

% Set range of kappa values of interest
idxrange = [24];

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
    for idx = idxrange
        for batchsize = batchrange
            for beta = momentumrange
                close all;
                
                fprintf("Start experiment for pfail = " + pfail / 10 + ", idx = " + ...
                    idx + ", batchsize = " + batchsize + ", momentum = " + beta + "\n");
                
                % loadzip(pfail / 10, idx);
                
                file = "zipcode_" + idx + "_pfail_" + pfail + ".mat";
                load(file);
                
                test_323_inner;
                
                fprintf("End experiment for pfail = " + pfail / 10 + ", idx = " + ...
                    idx + ", batchsize = " + batchsize + ", momentum = " + beta + "\n");
                
            end % End for
        end % End for
    end % End for
end % End for


