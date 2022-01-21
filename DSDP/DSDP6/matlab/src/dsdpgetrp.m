function [rp] = dsdpgetrp(A, X, b)
% Compute primal infeasibility AX - b

m = length(b);
rp = zeros(m, 1);

for i = 1:m
    
    rp(i) = trace(A{i} * X) - b(i);
    
end % End for


end % End function