function [step] = dsdpgetStepsize(S, dS, kappa, dkappa, tau, dtau, stepstrategy, alpha)

if stepstrategy == "linesearch"
    stepS = dsdpgetalpha(S, dS, alpha);
else
    stepS = dsdpgetalpha(S, dS);
end % End if

step = min([dtau / tau; dkappa / kappa]);

if step < 0
    step = abs(1 / step);
else
    step = 1e+03;
end % End if 

if stepS > 0
    step = min(step, stepS);
end % End if 

end % End function