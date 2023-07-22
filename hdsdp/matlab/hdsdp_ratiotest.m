function [alphamax] = hdsdp_ratiotest(S, dS)
% SDP ratio test

L = chol(S)';
B = L' \ dS;
B = L' \ B';
B = (B + B') / 2;

[lambda, delta, k] = hdsdp_lanczos(L', -dS);

alphamax = - 1 / eigs(B, 1, "smallestreal");

if isnan(alphamax)
    alphamax = - 1 / real(min(eig(full(B))));
end % End if 

end % End function