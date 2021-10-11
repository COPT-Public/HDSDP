% Outer test script for Zipcode Phase Retrieval Minibatch (3.2.2)
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
steprange = logspace(0, 2, 10);
nsteptotest = length(steprange);

% Set index of image 
digitimgs = [24];

% Set range of pfail values of interest (here we use pfail * 10 for 
% representation)
pfailrange = [2, 3];

% Set number of repeats
nTest = 1;

% Set optimization parameters
use_vm = 0;
beta = 0;
maxiter = 200;
show_info = false;

% Start experiment
% Set global random seed
rng(123);

for pfail = pfailrange
    for imgidx = digitimgs
        
        close all;

        fprintf("Start experiment for pfail = " + pfail / 10 + ", imgidx = " + ...
            imgidx + "\n");
        
        loadzip(pfail / 10, imgidx);
        
        file = "zipcode_" + imgidx + "_pfail_" + pfail + ".mat";
        load(file);

        test_322_inner;
        
        fprintf("Ended experiment for pfail = " + pfail / 10 + ", imgidx = " + ...
            imgidx + "\n");
    end % End for
end % End for


