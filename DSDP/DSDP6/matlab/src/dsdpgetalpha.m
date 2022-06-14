function [alphamax] = dsdpgetalpha(S, dS, alpha)
% The function computes the maximum stepsize to take such that
% S + alpha * dS is positive semi-definite 
if nargin == 3
    useline = true;
else
    useline = false;
end % End if

if useline
    alphamax = 1 / alpha;    
    while ~ dsdpIspsd(S + alphamax * dS) && alphamax > 1e-08
        alphamax = alphamax * 0.8;
    end % End if 
    if alphamax < 1e-08
        alphamax = 0.0;
    end % End if 
else
% Do Cholesky decomposition

R = dsdpperturbChol(S);
B = R' \ dS;
B = R' \ B';
B = (B + B') / 2;
% alphamax = - 1 / min(eig(B));
alphamax = - 1 / eigs(B, 1, "smallestreal");


% dsdplanczos(R, -dS);
% [eigmax, delta] = dsdplanczos(R, -dS);
% alphamax = 1 / (eigmax + delta);

if isnan(alphamax)
    alphamax = - 1 / real(min(eig(full(B))));
end % End if 
end % End if 

end % End function