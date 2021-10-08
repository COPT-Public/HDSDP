function [params] = setdefaultparams(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Inertial Model-based Stochastic Methods                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Oct 7th, 2020                                %
%                      Phase Retrieval Problem                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function provides a struct including default parameters for       %
% optimization                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                  %
%           A : problem data struct                                      %
% Output                                                                 %
%      params : a struct containing default parameters for optimization  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m, n] = size(data.A);
params.gamma = 0.1;
params.beta = 0.01;
params.init_x = randn(n, 1);
params.maxiter = 400;
params.tol = data.bestloss * 1.5;
params.early_stop = false;
params.batch = 1;
params.adp_param = 0;
params.alpha_0 = 1;
params.use_vm = 0.05;
params.show_info = true;
params.logscale = false;
params.optmethod = "ProxPt";

params.ninneriter = params.maxiter * m;
params.candidatemethod = ["SGD", "ProxLin", "ProxPt", "Vien"];

end % End function