function [obj] = getgoldobj(mu, csinv, b, asinv, M)
% Compute objective of golden linesearch
obj = mu * csinv + M * norm(b - mu * asinv, 1);
end % End function