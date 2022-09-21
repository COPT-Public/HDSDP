clear; clc; close all;

addpath ../src/;
addpath ../data/;

files = dir(fullfile("..", "data", "p_*"));
nfiles = length(files);
fnames = {files.name}';
maxmn = 1000;
minmn = 0;

diary("log_20220915.txt")
diary on
for i = 1:nfiles
    test_netlib(fullfile("..", "data", fnames{i}), 1000, maxmn, minmn);
end % End for
diary off;