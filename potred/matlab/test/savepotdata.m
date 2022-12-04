function [] = savepotdata(data, fname)

if nargin == 1
    fname = "pot-test";
end % End if

data = rmfield(data, 'objcon');
data = rmfield(data, 'lb');

save(fullfile('/Users/gaowenzhi/Desktop/gwz/potred/', 'src/potlp/potlp/', 'test/', fname) + ".mat", 'data');

end % End function