function [A, b, C, pscaler, dscaler] = dsdpPresolver(A, b, C)
% Implement the presolver (coefficient scaler) of DSDP

% Primal scaling
m = length(b);

pscaler = zeros(m, 1);
nnormarray = zeros(m + 1, 1);
nnormarray(m + 1) = norm(C, 'fro');
for j = 1:m
    nnormarray(j) = norm(A{j}, 'fro');
end % End for

for i = 1:m
    pscaler(i) = sqrt(nnormarray(i) * abs(b(i)));

    if pscaler(i) == 0
        pscaler(i) = 1;
    end % End if 
    
    A{i} = A{i} / pscaler(i);
    b(i) = b(i) / pscaler(i);
    % nnormarray(i) = nnormarray(i) / pscaler(i);
end % End for

% nnormarray = nnormarray(nnormarray > 0);
% Only one block 
dscaler = 1.0;
% dscaler = sqrt(max(nnormarray) * min(nnormarray));

% for i = 1:m
%     A{i} = A{i} / dscaler;
% end % End for 
% 
% C = C / dscaler;

end % End function