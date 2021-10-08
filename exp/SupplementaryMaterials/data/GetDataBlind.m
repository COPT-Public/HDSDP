function [U, V, b] = GetDataBlind(x, m, kappa, pfail)

% Generate data
[U, ~] = GetData(x, m, kappa, pfail);
[V, ~] = GetData(x, m, kappa, pfail);

% Compute b
b = (U * x) .* (V * x) + (rand(m, 1) < pfail).* randn(m, 1) * 5;

end % End function