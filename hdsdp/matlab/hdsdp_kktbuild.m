function [M, asinv, asinvrdsinv, asinvcsinv, csinv, csinvrdcsinv] =...
    hdsdp_kktbuild(As, C, S, Rd)
% Build the KKT system
m = length(As);
Sinv = inv(S);

M = zeros(m, m);
asinv = zeros(m, 1);
asinvrdsinv = zeros(m, 1);

csinv = trace(C * Sinv); %#ok
csinvrdcsinv = Rd * trace(Sinv * C * Sinv); %#ok

for i = 1:m
    asinv(i) = trace(As{i} * Sinv); %#ok
    asinvcsinv(i) = trace(As{i} * Sinv * C * Sinv); %#ok
    asinvrdsinv(i) = Rd * trace(Sinv * As{i} * Sinv); %#ok
    for j = 1:i
        M(i, j) = trace(As{i} * Sinv * As{j} * Sinv); %#ok
        M(j, i) = trace(As{i} * Sinv * As{j} * Sinv); %#ok
    end % End for
end % End for

end % End function