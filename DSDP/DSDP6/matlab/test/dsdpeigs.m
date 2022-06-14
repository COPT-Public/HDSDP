clear;
clc;

load Aeig.mat;
Annz = (A ~= 0);

Anzrow = sum(Annz, 1);
