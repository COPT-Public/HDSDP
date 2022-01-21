function [step] = dsdpgetStepsize(S, dS, kappa, dkappa, tau, dtau, stepstrategy, alpha)

if stepstrategy == "linesearch"
    stepS = dsdpgetalpha(S, dS, alpha);
else
    stepS = dsdpgetalpha(S, dS);
end % End if

step = min([dtau / tau; dkappa / kappa]);

if step < 0
    step = min(1, abs(alpha / step));
else
    step = 1.0;
end % End if 

if stepS > 0
    step = min([step, alpha * stepS, 1]);
end % End if 

end % End function