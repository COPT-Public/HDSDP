function [ATy] = hdsdp_aty(A, y)

ATy = A{1} * y(1);
for i = 2:length(y)
    ATy = ATy + A{i} * y(i);
end % End for

end % End function