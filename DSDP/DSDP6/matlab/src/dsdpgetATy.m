function [ATy] = dsdpgetATy(A, y)
% This function computes ATy = A{1} * y(m)+ ... + A{m} * y(m)

% Get size
m = length(y);
[n, ~] = size(A{1});

% Prepare solution
ATy = sparse(n, n);

% Assemble M
for i = 1:m
   
    ATy = ATy + y(i) * A{i};
    
end % End parfor

% End function