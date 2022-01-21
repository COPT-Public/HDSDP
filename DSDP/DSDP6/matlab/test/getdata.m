clear;


% diary infeas_bench.txt
% diary on;

problems = [ "control3.dat-s", "gpp124-2.dat-s", "hinf1.dat-s", "hinf6.dat-s",...
    "maxG55.dat-s", "mcp500-1.dat-s", "qpG51.dat-s", "truss2.dat-s", ...
    "control4.dat-s", "gpp124-3.dat-s", "hinf10.dat-s", "hinf7.dat-s",...
    "maxG60.dat-s", "mcp500-2.dat-s", "ss30.dat-s", "truss3.dat-s", ...
    "control5.dat-s", "gpp124-4.dat-s", "hinf11.dat-s", "hinf8.dat-s",...
    "mcp100.dat-s", "mcp500-3.dat-s", "theta1.dat-s", "truss4.dat-s",...
    "arch0.dat-s", "control6.dat-s", "gpp250-1.dat-s", "hinf12.dat-s",...
    "hinf9.dat-s", "mcp124-1.dat-s", "mcp500-4.dat-s", "theta2.dat-s",...
    "truss5.dat-s", "arch2.dat-s", "control7.dat-s", "gpp250-2.dat-s",...
    "hinf13.dat-s", "infd1.dat-s", "mcp124-2.dat-s", "qap10.dat-s",...
    "theta3.dat-s", "truss6.dat-s", "arch4.dat-s", "control8.dat-s",...
    "gpp250-3.dat-s", "hinf14.dat-s", "infd2.dat-s", "mcp124-3.dat-s",...
    "qap5.dat-s", "theta4.dat-s", "truss7.dat-s", "arch8.dat-s",...
    "control9.dat-s", "gpp250-4.dat-s", "hinf15.dat-s", "infp1.dat-s",...
    "mcp124-4.dat-s", "qap6.dat-s", "theta5.dat-s", "truss8.dat-s",...
    "control1.dat-s", "equalG11.dat-s", "gpp500-1.dat-s", "hinf2.dat-s", ...
    "infp2.dat-s", "mcp250-1.dat-s", "qap7.dat-s", "theta6.dat-s", ...
    "control10.dat-s", "equalG51.dat-s", "gpp500-2.dat-s", "hinf3.dat-s", ...
    "maxG11.dat-s", "mcp250-2.dat-s", "qap8.dat-s", "thetaG11.dat-s", ...
    "control11.dat-s", "gpp100.dat-s", "gpp500-3.dat-s", "hinf4.dat-s", ...
    "maxG32.dat-s", "mcp250-3.dat-s", "qap9.dat-s", "thetaG51.dat-s", ...
    "control2.dat-s", "gpp124-1.dat-s", "gpp500-4.dat-s", "hinf5.dat-s", ...
    "maxG51.dat-s", "mcp250-4.dat-s", "qpG11.dat-s", "truss1.dat-s"];

ntest = length(problems);
iscorrect = zeros(ntest, 1);

for k = 1:ntest
    
    % [At,b,c,K] = readsdpa(fullfile('../benchmark/sdplib/', 'gpp100.dat-s'));
    [At,b,c,K] = readsdpa(fullfile('../benchmark/sdplib/', problems(k)));
    % [At,b,c,K] = readsdpa(fullfile('../benchmark/DIMACS/', 'prob_1_2_0.dat-s'));
    % load(fullfile("../benchmark", "bm1.mat"));
    % load(fullfile("../benchmark", "infeas", "infeas_messy_20_10", "infeas_messy_20_10_" + k + ".mat"));
    
    % fprintf("\n Test %s \n", problems(k));
    try
        [nsqr, ~] = size(At);
        A = At';
    catch
        [~, nsqr] = size(A);
    end % End try
    
    if nsqr > 100000
        continue;
    end % End if
    
    A = A(:, 1:nsqr);
    [m, n] = size(A);
    p = m;
    n = sqrt(n);
    Amat = cell(p, 1);
    
    try
        for i = 1:p
            Amat{i} = reshape(A(i, :), n, n);
            Amat{i} = (Amat{i} + Amat{i}') / 2;
        end % End for
    catch
        continue;
    end % End try
    C = reshape(c(1:nsqr), n, n);
    
    fprintf("Done with %s \n", problems(k));
    eval('save ' + problems(k) + '.mat  Amat b C;');
    
end % End for

