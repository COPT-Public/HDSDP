function [isray] = dsdpcheckDualRay(A, b, dy)
% Check if a direction implies a dual improving ray

isray = false;
if b' * dy < 0
    return;
end % End if

if dsdpIspsd(dsdpgetATy(A, -dy))
    isray = true;
end % End if 

end % End function