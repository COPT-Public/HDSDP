function [merit] = dsdpgetMeritValueRlx(b, y, mu, L, bound)
% Get the value of the merit function

merit = b' * y + mu * 2 * sum(log(diag(L)));
merit = merit + mu * (sum(log(bound - y)) + sum(log(bound + y)));
merit = - merit;

end % End function