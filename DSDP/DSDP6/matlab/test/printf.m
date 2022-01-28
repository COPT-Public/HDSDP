function [] = printf(num)

if issparse(num)
    num = full(num);
end % End if 

fprintf("    %30.30e \n", num);


end % End function