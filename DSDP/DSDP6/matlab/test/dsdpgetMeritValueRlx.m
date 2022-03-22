function [merit] = dsdpgetMeritValueRlx(b, y, mu, L)
% Get the value of the merit function

merit = b' * y + mu * 2 * sum(log(diag(L)));
merit = merit + mu * (sum(log(1e+07 - y)) + sum(log(1e+07 + y)));
merit = - merit;

end % End function