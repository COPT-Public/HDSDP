% Test script for Zipcode Phase Retrieval Image Recovery (3.2.4)
clear;
clc;
close all;

addpath(fullfile("..", "opt"));
addpath(fullfile("..", "data"));

idx = 2;
pfail = 3;
loadzip(pfail / 10, idx);
file = "zipcode_" + idx + "_pfail_" + pfail + ".mat";
load(file);

bestloss = data.bestloss;
params = setdefaultparams(data);
params.use_vm = 0;
params.maxiter = 200;
params.show_info = true;
params.early_stop = false;
params.beta = 0;
params.gamma = 1/sqrt(params.maxiter * 768 / params.batch);
% params.gamma = 1;
params.adp_param = 0;
params.alpha_0 = 100;
params.logscale = true;
params.show_info = false;

nTest = 1;
batchrange = [4, 8, 16, 32, 48, 64];

params.init_x = 2 * randn(256, 1) + data.optx;

for i = 1:nTest
    
    params.gamma = 1/sqrt(params.maxiter * 768 / params.batch);
    
%     % Test SGD
%     params.optmethod = "SGD";
%     [sgdsol, sgdinfo] = proxopt(data, params);
    
    % Test Proximal linear
    params.optmethod = "ProxLin";
    [proxlinsol, proxlininfo] = proxopt(data, params);
    
    
end % End for

for i = 1:11
    subplot(7, 11, i);
    imshow(reshape(proxlininfo.imgres(i, :), 16, 16)')
end % End for


for i = 1:length(batchrange)
    
    params.batch = batchrange(i);
    % params.alpha_0 = sqrt(params.batch);
    params.gamma = 1/sqrt(params.maxiter * 768 / params.batch);
    
    for j = 1:nTest
        
%         % Test SGD
%         params.optmethod = "SGD";
%         [sgdsolb, sgdinfob] = proxopt(data, params);
        
        % Test Proximal linear
        params.optmethod = "ProxLin";
        [proxlinsolb, proxlininfob] = proxopt(data, params);
        
    end % End for
    
    for k = 1:11
        subplot(7, 11, k + 11 * i);
        imshow(reshape(proxlininfob.imgres(k, :), 16, 16)')
    end % End for
    
    fprintf("Done for batchsize " + i + " \n");
    
end % End for

save("zipcode_" + idx + "_pfail_" + pfail + "_image.mat");
savefig("zipcode_" + idx + "_pfail_" + pfail + "_image");
