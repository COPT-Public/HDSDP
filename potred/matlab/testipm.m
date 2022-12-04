clear; clc; close all;

fname = fullfile('data', 'p_tuff.mps');

data = preprocess(fname);

A = data.A;
b = data.b;
% b = b / (norm(b, 1) + 1);
c = data.c;
% c = c / (norm(c, 1) + 1);

tic;
[x, y, s, kappa, tau] = hsdipm(A, b, c);
toc
% tic;
% [x, y, s, kappa, tau] = hsdipm2(A, b, c);
% toc