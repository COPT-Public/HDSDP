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
ntest = 1; % length(problems);
iscorrect = zeros(ntest, 1);

for k = 1:ntest
    
    [At,b,c,K] = readsdpa(fullfile('../benchmark/sdplib/', 'theta1.dat-s'));
    % [At,b,c,K] = readsdpa(fullfile('../benchmark/sdplib/', problems(k) + '.dat-s'));
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
    
    % nsqr = 10000;
    
    A = A(:, 1:nsqr);
    [m, n] = size(A);
    p = min(m, m);
    n = sqrt(n);
    Amat = cell(p, 1);
    
    for i = 1:p
        Amat{i} = reshape(A(i, :), n, n);
        Amat{i} = (Amat{i} + Amat{i}') / 2;
    end % End for
    
    dsdpParam = dsdpgetParam();
    
    objacc = zeros(ntest, 1);
    solacc = zeros(ntest, 1);
    niter = zeros(ntest, 1);
    
    
    %     A2 = zeros(p, n);
    %     for i = 1:p
    %         A2(i, :) = sqrt(diag(Amat{i}));
    %     end % End for
    
    C = reshape(c(1:nsqr), n, n);
    % C = C / norm(C, 1);
    b = b(1:p);
    warning off;
    %     cvx_begin sdp
    %     cvx_solver mosek
    %     variable X(n, n) symmetric
    %     dual variable S1
    %     dual variable y1
    %     minimize trace(C * X)
    %
    %     trace(Amat{1} * X) == b(1) : y1;
    %     for i = 2:p
    %         trace(Amat{i} * X) == b(i);
    %     end % End for
    %
    %     X >= 0 : S1;
    %     cvx_end
    
    pinfeas = zeros(p, 1);
    
    % if cvx_status == "Solved"
    try
        [X2, S2, y, kappa, tau, reason] = dSDPKappaTauPds(Amat, b, C, dsdpParam);
        
        fprintf("\n\nSummary \n");
        %         fprintf("Objective Error: %3.e\n", abs(b' * y - trace(C * X)));
        %         fprintf("Solution  Error: %3.e\n", norm(S1 - S2, "fro"));
        %
        %         objacc(k) = abs(b' * y - trace(C * X)) / trace(C * X);
        %         solacc(k) = norm(S1 - S2, "fro") / norm(S1, "fro");
        
        for i = 1:p
            pinfeas(i) = trace(Amat{i} * X2) - b(i);
        end % End for
        
        Rd = dsdpgetATy(Amat, y) + S2 - C;
        
        bn1 = norm(b, 1);
        cn1 = sum(sum(abs(C)));
        err1 = norm(pinfeas) / (1 + bn1);
        err2 = max(0, - eigs(X2, 1, 'smallestreal')) / (1 + bn1);
        err3 = norm(Rd, 'fro') / (1 + cn1);
        err4 = max(0, - eigs(S2, 1, 'smallestreal')) / (1 + cn1);
        err5 = (trace(C * X2) - b' * y) / (1 + abs(trace(C * X2)) + abs(b' * y));
        err6 = trace(X2 * S2) / (1 + abs(trace(C * X2)) + abs(b' * y));
        fprintf("DIMACS Error: %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e \n", err1, err2, err3, err4, err5, err6);
        
        if reason == "DSDP_PRIMAL_INFEASIBLE_DUAL_UNBOUNDED"
            iscorrect(k) = 1;
        else
            1;
        end % End if
        
        
    catch
        iscorrect(k) = 0;
    end % End try
    
    % end % End if
end % End for

diary off;
