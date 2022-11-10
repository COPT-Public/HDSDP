clear; clc; close all;

fname = fullfile('data', 'p_SCSD1.SIF.mps');

data = preprocess(fname);

A = data.A;
b = data.b;
% b = b / (norm(b, 1) + 1);
c = data.c;
% c = c / (norm(c, 1) + 1);

[x, y, s, kappa, tau] = hsdipm(A, b, c);