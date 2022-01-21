function [merit] = dsdpgetMeritValue(b, y, mu, L)
% Get the value of the merit function

merit = b' * y + mu * 2 * sum(log(diag(L)));

end % End function