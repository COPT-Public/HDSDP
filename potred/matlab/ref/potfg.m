function [phi, gphi] = potfg(rho, f, g, m, x, s, kappa, tau)

phi = rho * log(f) - sum(log(x)) - sum(log(s)) - log(kappa) - log(tau); 
g2 = [zeros(m, 1); 1./x; 1./s; 1/kappa; 1/tau];
gphi = rho * g / f - g2;    
    
end % 