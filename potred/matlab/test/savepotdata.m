function [] = savepotdata(data, fname)

if nargin == 1
    fname = "pot-test";
end % End if

data = rmfield(data, 'objcon');
data = rmfield(data, 'lb');

save(fname + ".mat", 'data');

end % End function