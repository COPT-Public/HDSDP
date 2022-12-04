clear; clc; close all;

cd /Users/gaowenzhi/Desktop/gwz/potred/matlab/test;
addpath ../src/;
addpath ../data/;

files = dir(fullfile("..", "data", "p_*"));
nfiles = length(files);
fnames = {files.name}';
maxmn = 1e+10;
minmn = 0;

if nfiles <= 1
    return;
end % End if

diary off;

log_path = "log_20221116";
mkdir(log_path);

fileID = fopen('all_1116.txt', 'w');
fprintf(fileID, "| %30s |   pObj   |   dObj   | pInfeas | dInfeas | relGap  |  Time | Status \n",...
        "Instance");

for i = 1:nfiles
    n = fnames{i};
    dname = fullfile(log_path, n(3:end-8) + ".txt");
%     delete(dname);
    diary(dname);
    diary on
    try 
        data = preprocess(fullfile("..", "data", fnames{i}));
%         savepotdata(data, n(3:end-4));
        test_netlib(fullfile("..", "data", fnames{i}), 1000000, 600.0, maxmn, minmn, fileID);
    catch
        
    end % End try
    diary off;
end % End for
