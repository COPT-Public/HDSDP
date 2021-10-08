function [sol, info] = proxopt(data, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Inertial Model-based Stochastic Methods.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Oct 7th, 2020                                %
%                      Phase Retrieval Problem                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function integrates three model-based methods to solve            %
% Phase Retrieval Problem                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                                 %
% data : a struct containing data for optimization                       %
%          A, b : data in optimization                                   %
% params : parameters for the optimization algorithm                     %
%   gamma, beta : optimization parameters, gamma for stepsize parameter  %
%                 and beta is for momentum parameter                     %
%        init_x : initial starting point of algorithm                    %
%       maxiter : maximum number of iterations (epochs) allowed          %
%           tol : tolerance allowed to end the algorithm prematurely     %
%    early_stop : whether to stop optimization when tolerance is reached %
%     adp_param : parameters for adaptive method                         %
%        use_vm : parameters for variable metric                         %
%     show_info : whether to display optimization progress               %
%     optmethod : which model-based algorithm to use                     %
% Output:                                                                %
%           sol : an array storing the last/best solution after the      %
%                 optimization procedure                                 %
%          info : an array storing the information related to the        %
%                 optimization procedure                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack data arrays
A = data.A;
b = data.b;

% Unpack parameter arrays
gamma = params.gamma;
beta = params.beta;
init_x = params.init_x;
maxiter = params.maxiter;
tol = params.tol;
early_stop = params.early_stop;
batch = params.batch;
adp_param = params.adp_param;
alpha_0 = params.alpha_0;
use_vm = params.use_vm;
show_info = params.show_info;
optmethod = params.optmethod;

% Check data validity
if show_info
    disp("Model-based optimization routine start on " + date);
    disp("Checking data validity");
end % End if

[m, n] = size(A);

if m == 0 || n == 0
    error("Empty data array");
end % End if

if ~ (size(b) == m)
    error("Sizes of A and b mismatch");
end % End if

if min(b) < 0
    error("b must be non-negative");
end % End

if show_info
    disp("Data validity is successfully checked");
end

% Check parameter validity
if show_info
    disp("Checking parameter validity");
end % End if

if gamma < 0 || beta < 0
    error("Optimization parameter gamma and beta must be non-negative");
end % End if

if maxiter <= 0
    error("Number of epochs shold at least be 1");
end % End if

if batch > m
    error("Batchsize exceeded: please use batchsize that " +...
        "is less than " + m);
end % End if

if show_info
    disp("Parameters successfully checked");
end % End if

% Check solution method
solmethod = "Naive";
if adp_param > 0
    solmethod = "Adaptive";
end % End if

if use_vm
    if solmethod == "Adaptive"
        warning("Variable metric will turn off adaptive method");
        adp_param = 0;
    end % End if
    solmethod = "Variable Metric";
end % End if

% Start the optimization
if show_info
    format short;
    disp("Optimization information summary:");
    disp("Data size: " + m + " samples of dimension " + n);
    disp("Optimization method: " + optmethod);
    disp("Solution method: " + solmethod);
    disp("Maximum number of iterations: " + maxiter);
    
    if early_stop
        disp("Early Stopping is turned on with tolerance " + tol);
    end % End if
    
    obj = sum(abs((A * init_x).^2 - b)) / m;
    disp("***********************************************************");
    disp("****               Optimization Start                  ****");
    disp("***********************************************************");
    disp("- Initial - Obj: " + obj + " - Best obj: " + obj);
    
end % End if

if optmethod == "SGD"
    [sol, info] = proxsgd(A, b, 1 / gamma, beta, init_x, maxiter, tol,...
        early_stop, batch, adp_param, alpha_0, use_vm, show_info);
elseif optmethod == "ProxLin"
    if batch == 1
        [sol, info] = proxlin(A, b, 1 / gamma, beta, init_x, maxiter, tol,...
            early_stop, adp_param, alpha_0, use_vm, show_info);
    else
        [sol, info] = proxlinbatch(A, b, 1 / gamma, beta, init_x, maxiter, tol,...
            early_stop, batch, adp_param, alpha_0, show_info);
    end % End if
else
    if batch == 1
        [sol, info] = proxpt(A, b, 1 / gamma, beta, init_x, maxiter, tol,...
            early_stop, adp_param, alpha_0, use_vm, show_info);
    else
        [sol, info] = proxptbatch(A, b, 1 / gamma, init_x, maxiter, tol,...
            early_stop, batch, adp_param, alpha_0, show_info);
    end % End if
end % End if

% Collect information in the info struct
info.solmethod = solmethod;
info.optmethod = optmethod;

if params.logscale
    info.bestobjs = log10(info.bestobjs);
end % End if

if show_info
    disp("***********************************************************");
    disp("****               Optimization Ends                   ****");
    disp("***********************************************************");
end % End

end % End function