clear; clc; close all;

addpath ../src/;
addpath ../data/;

files = dir(fullfile("..", "data", "p_*"));
nfiles = length(files);
fnames = {files.name}';
maxmn = 500;
minmn = 0;

diary off;

for i = 1:nfiles
    n = fnames{i};
    diary("log_20221027" + n(3:end-8) + ".txt");
    diary on
    test_netlib(fullfile("..", "data", fnames{i}), 500000, maxmn, minmn);
    diary off;
end % End for
