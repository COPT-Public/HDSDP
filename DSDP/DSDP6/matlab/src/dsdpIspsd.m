function [retcode] = dsdpIspsd(A)

retcode = true;

% if eigs(A, 1, 'smallestreal') < 0
%     retcode = false;
% end % End if 

try
    lchol(A);
catch
    retcode = false;
end % End try

end % End function