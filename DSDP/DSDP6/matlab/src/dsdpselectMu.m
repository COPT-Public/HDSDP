function [newmu] = dsdpselectMu(A, S, muk, dy, dy1, backwardnewton, ismukfeas)
% Take primal step and find new mu parameter
% dy = dy1 / muk + dy2

if ismukfeas
    dnewton = dsdpgetATy(A, dy1) / muk;
    alphamu = dsdpgetalpha(backwardnewton, dnewton);
    if alphamu < 0
        alphamu = 1.0;
    end % End if
    newmu = muk / (1 + 0.95 * alphamu);
else
    dS = dsdpgetATy(A, dy);
    alphap = dsdpgetalpha(S, dS);
    assert(alphap > 0);
    Shat = S + 0.95 * alphap * dS;
    % assert(dsdpIspsd(Shat));
    dS = - alphap * dsdpgetATy(A, dy1) / muk;
    alphamu = dsdpgetalpha(Shat, dS, 1.0);
    if alphamu < 0
        alphamu = 1.0;
    end % End if 
    newmu = (alphap * muk) / (1 + alphamu) + (1 - alphap) * muk;
end % End if 

% End function