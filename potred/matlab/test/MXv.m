function [mxv] = MXv(v, coneidx, AT, A, ATA, x, f, g, rho, scale)
% Compute the Hessian-vector product for negative curvature computation

n = length(v);
ncone = length(coneidx);

if scale
    nrm = norm(x(coneidx));
    xp = x(coneidx) / nrm;
    v(coneidx) = v(coneidx) - (xp' * v(coneidx)) * xp;
    % Set up u1
    u1 = zeros(n, 1); u1(coneidx) = v(coneidx);
    
    % Set up u2
    u2 = v; 
    u2(coneidx) = u2(coneidx) .* x(coneidx);
    u2 = AT * (A * u2);
    u2(coneidx) = u2(coneidx) .* x(coneidx);
    u2(coneidx) = u2(coneidx) - (xp' * u2(coneidx)) * xp;
    
    % Set up u3
    g(coneidx) = g(coneidx) .* x(coneidx);
    u3 = (g' * v) * g;
    u3(coneidx) = u3(coneidx) - (xp' * u3(coneidx)) * xp;
    mxv = f * u1 + rho * u2 - rho / f * u3;
else
    % Not recommended. For testing purpose
    xcone = x(coneidx);
    v(coneidx) = v(coneidx) - sum(v(coneidx)) / ncone;
    u1 = zeros(n, 1); u1(coneidx) = v(coneidx);
    u1(coneidx) = u1(coneidx) ./ (xcone.^2);
    u1(coneidx) = u1(coneidx) - sum(u1(coneidx)) / ncone;
    u2 = AT * (A * v);
    u2(coneidx) = u2(coneidx) - sum(u2(coneidx)) / ncone;
    u3 = (g' * v) * g;
    u3(coneidx) = u3(coneidx) - sum(u3(coneidx)) / ncone;
    mxv = f^2 * u1 + rho * f * u2 - rho * u3;
end % End if

end % End function