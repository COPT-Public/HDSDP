clear;
clc;

m = 3;
n = 1000;
ntest = 1;
sps = 0.001;

rng(24);

dsdpParam = dsdpgetParam();

objacc = zeros(ntest, 1);
solacc = zeros(ntest, 1);
niter  = zeros(ntest, 1);

for k = 1:ntest
    A = cell(m, 1);
    
    for i = 1:m
        A{i} = sprandsym(n, sps);
    end % End for
    
    b = rand(m, 1);
    C = sparse(1:n, 1:n, rand(n, 1));
    % C = C' * C;
    
%     cvx_begin sdp
%     cvx_solver mosek
%     variable X(n, n) symmetric
%     dual variable S1
%     dual variable y1
%     minimize trace(C * X)
%     
%     trace(A{1} * X) == b(1) : y1;
%     for i = 2:m
%         trace(A{i} * X) == b(i);
%     end % End for
%     
%     X >= 0 : S1;
%     cvx_end
    
    pinfeas = zeros(m, 1);
   
    % if cvx_status == "Solved"
    tic;
        [X2, S2, y, kappa, tau, iter] = dSDPKappaTauPds(A, b, C, dsdpParam);
    toc
        fprintf("\n\nSummary \n");
        niter(k) = iter;
        
        for i = 1:m
            pinfeas(i) = trace(A{i} * X2) - b(i);
        end % End for 
        
        Rd = dsdpgetATy(A, y) + S2 - C;
        
        bn1 = norm(b, 1);
        cn1 = sum(sum(abs(C)));
        err1 = norm(pinfeas) / (1 + bn1);
        err2 = max(0, - eigs(X2, 1, 'smallestreal')) / (1 + bn1);
        err3 = norm(Rd, 'fro') / (1 + cn1);
        err4 = max(0, - eigs(S2, 1, 'smallestreal')) / (1 + cn1);
        err5 = (trace(C * X2) - b' * y) / (1 + abs(trace(C * X2)) + abs(b' * y));
        err6 = trace(X2 * S2) / (1 + abs(trace(C * X2)) + abs(b' * y));
        fprintf("DIMACS Error: %e %e %e %e %e %e \n", err1, err2, err3, err4, err5, err6);
        
    % end % End if
end % End for
% 
% fprintf("ObjAcc: %5.3e \n", mean(objacc));
% fprintf("SolAcc: %5.3e \n", mean(solacc));
% fprintf("nIter: %5.3e \n", mean(niter));


