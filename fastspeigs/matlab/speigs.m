%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPEIGS: Sparse eigen-value decomposition v1.0.0, Aug 27th, 2022
% Wenzhi Gao, Shanghai University of Finance and Economics
% 
% The routine implements utility that performs eigen-decomposition on
% extremely sparse matrices that usually arise from semi-definite
% programming problems. Specially structures are exploited to reduce the
% computation time.
% 
% Please refer to the speigs/doc for detailed documents
% 
% Usage:
% 
%       [V, e] = speigs(A, opts);
%
%   A    : n * n sparse matrix to be decomposed, ONLY lower triangular is
%          needed
%   opts : options for the algorithm
%   V    : n * n matrix containing eigen-vectors
%   e    : n * 1 array containing eigen-values
% 
% Parameters:
%   opts.gthresh : When the submatrix size exceeds gthresh * n, it is
%                  considered as a general matrix 
%   opts.tol     : Eigen-values whose abs <= tol will be considered 0
%   opts.quiet   : If quiet, no information will be printed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
