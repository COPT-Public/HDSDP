function [X, y] = dsdpPostsolver(X, y, pscaler, dscaler)
% Post-solver for DSDP

X = X * dscaler;
y = y .* pscaler;

end % End function
