clear; clc; close all;

files = dir(fullfile("/Users/gaowenzhi/Desktop/potential-reduction/src/potlp/potlp/test/stform", "*.mps"));
nfiles = length(files);
fnames = {files.name}';

% model = gurobi_read(fullfile('/Users/gaowenzhi/Desktop/potential-reduction/src/potlp/potlp/test/stform', 's_ken-07.mps'));
% A = model.A;
% b = model.rhs;
% c = model.obj;
% 
% [D, E, A] = l2scale(A);
% b = b ./ D;
% c = c ./ E;
% [x, y, s] = dualpot(A, b, c, 1);
% 
% return;

for i = 1:nfiles
    
    n = fnames{i};
    
    model = gurobi_read(fullfile('/Users/gaowenzhi/Desktop/potential-reduction/src/potlp/potlp/test/stform', fnames{i}));
    A = model.A;
    b = model.rhs;
    c = model.obj;
    
    [D, E, A] = l2scale(A);
    b = b ./ D;
    c = c ./ E;
    
    model.A = A;
    model.rhs = b;
    model.obj = c;
    
    [x, y, s] = dualpot(A, b, c);

    try
        [x, y, s] = dualpot(A, b, c);
        param.OutputFlag = 0;
        sol = gurobi(model, param);
        fprintf("%20s %5e \n", fnames{i}, abs(c' * sol.x - b' * y) / (abs(sol.objval) + 1));
    catch
        fprintf("%20s %5e \n", fnames{i}, "-1");
    end % End try

end % End for