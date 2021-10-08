function [A, b] = loadzip(pfail, idx)
% This function will return a phase retrieval dataset based on mnist

% Phase Retrieval Dataset Preparation
load zipcode.mat;

% Get one image
img = data(idx, 1:256)';

% Generate Hadamard matrix
Hmat = hadamard(256) / 16;

% Generate A and b
A = [Hmat * diag(2 * (randn(256, 1) >= 0) - 1); 
    Hmat * diag(2 * (randn(256, 1) >= 0) - 1); 
    Hmat * diag(2 * (randn(256, 1) >= 0) - 1)]; 
b = ((A * img).^2) .* (rand(768, 1) >= pfail);

clear data;

% Save data
data.A = A;
data.b = b;
data.optx = img;
data.bestloss = sum(abs((A * img).^2 - b)) / 768;
save("zipcode_" + idx + "_pfail_" + 10 * pfail + ".mat", "data");

end % End for


