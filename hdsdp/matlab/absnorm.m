function [onenrm] = absnorm(A)

onenrm = full(sum(sum(abs(A))));

end % End function