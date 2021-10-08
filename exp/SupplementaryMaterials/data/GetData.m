function [A, b] = GetData(x, m, kappa, pfail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inertial Model-based Stochastic Methods for    %
%       Nonsmooth Nonconvex Optimization           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates dataset for phase retrieval
% Input: 
%        x: real model parameter of dimension n
%        m: number of samples generated
%
% Output:
%        A: an m * n matrix of m samples
%        b: b_i is computed by b_i = (a_i' * x)^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get dimension of x
n = length(x);

% Generate A = Q * D
Q = randn(m, n);
D = diag(linspace(1/kappa, 1, n));
A = Q * D;

% Compute b
b = (A * x).^2 + (rand(m, 1) < pfail).* randn(m, 1) * 5;

end % End function



