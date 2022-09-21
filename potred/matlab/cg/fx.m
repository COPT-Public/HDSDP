function [fval] = fx(A, b, x)

fval = 0.5 * x' * A * x - b' * x;

end % End function