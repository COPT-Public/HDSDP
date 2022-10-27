clear; clc; close all;

addpath ../src/;
addpath ../data/;

files = dir(fullfile("..", "data", "p_*"));
nfiles = length(files);
fnames = {files.name}';
maxmn = 400;
minmn = 0;

diary off;

log_path = "log_20221027";
mkdir(log_path);

for i = 1:nfiles
    n = fnames{i};
    dname = fullfile(log_path, n(3:end-8) + ".txt");
%     delete(dname);
    diary(dname);
    diary on
    try 
        test_netlib(fullfile("..", "data", fnames{i}), 150000, 60.0, maxmn, minmn);
    catch
        
    end % End try
    diary off;
end % End for
