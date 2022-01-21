function [pinv] = dsdpInvPerm(p)

% Get the inverse of some permutation
n = length(p);
pinv = 1:n;
pinv(p) = 1:n;

end % End function