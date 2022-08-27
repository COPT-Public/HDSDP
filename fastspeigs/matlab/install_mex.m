clear; clc; close all;

version = "1.0.0";
fprintf("Installing SPEIGS %s \n", version);

addpath(".");
copyfile(fullfile("build", "speigs.*"), ".");
test_speigs;

fprintf("\nInstallation is done \n");
