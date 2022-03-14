function [newmu] = dsdpselectMuRlx(A, sl, su, S, muk, dy, dy1, backwardnewton, backwardnewtonlb, backwardnewtonub, ismukfeas, gap)
% Take primal step and find new mu parameter
% dy = dy1 / muk + dy2

if ismukfeas
    dnewton = dsdpgetATy(A, dy1) / muk;
    dnewtonub = dy1 / muk;
    dnewtonlb = -dnewtonub;
    
    alphamu = dsdpgetalpha(backwardnewton, dnewton);
    if alphamu < 0
        alphamu = 1.0;
    end % End if
    
    alphamu2 = -1 / min(dnewtonlb ./ backwardnewtonlb);
    if alphamu2 < 0
        alphamu2 = 1.0;
    end % End if 
    
    alphamu3 = -1 / min(dnewtonub ./ backwardnewtonub);
    if alphamu3 < 0
        alphamu3 = 1.0;
    end % End if 
    
    alphamu = min([alphamu, alphamu2, alphamu3]);
    newmu = muk / (1 + 0.95 * alphamu);
    
else
    dS = dsdpgetATy(A, dy);
    dls = - dy; 
    dus = dy;
    alphap = dsdpgetalpha(S, dS);
    alphap2 = -1 / min(dls ./ sl);
    alphap3 = -1 / min(dus ./ su);
    
    if alphap < 0
        alphap = 1.0;
    end % End if 
    
    if alphap2 < 0
        alphap2 = 1.0;
    end % End if 
    
    if alphap3 < 0
        alphap3 = 1.0;
    end % End if 
    
    alphap = min([alphap, alphap2, alphap3]);
    
    assert(alphap > 0);
    Shat = S + 0.95 * alphap * dS;
    slhat = sl + 0.95 * alphap * dls;
    suhat = su + 0.95 * alphap * dus;
    % assert(dsdpIspsd(Shat));
    dS = - alphap * dsdpgetATy(A, dy1) / muk;
    dls = - alphap * dy1 / muk;
    dus = alphap * dy1 / muk;
    alphamu = dsdpgetalpha(Shat, dS);
    alphamu1 = -1 / min(dls ./ slhat);
    alphamu2 = -1 / min(dus ./ suhat);
    
    if alphamu < 0
        alphamu = 1.0;
    end % End if 
    
    if alphamu1 < 0
        alphamu1 = 1.0;
    end % End if 
    
    if alphamu2 < 0
        alphamu2 = 1.0;
    end % End if
    
    alphamu = min([alphamu, alphamu1, alphamu2]);
    newmu = (alphap * muk) / (1 + alphamu) + (1 - alphap) * gap;
    
end % End if 

end % End function