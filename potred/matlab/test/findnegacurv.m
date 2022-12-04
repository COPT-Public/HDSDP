function [v, vbef] = findnegacurv(x, m, coneidx, projidx, rho, g, f, ATA, AT, A, vstart, method)
% Find negative curvature of the Hessian matrix
% Hess =  rho * (- (g * g') / f + ATA) + diag(f * d)
% over the subspace e' * v = 0

n = length(x);
ncone = length(coneidx);
nproj = length(projidx);
e = ones(ncone, 1);
e(projidx - min(coneidx) + 1) = 0;
e = 1 - e;

% Filter the non-basic variables
% bid = find(x(coneidx) > 5e-04);
% [~, bid] = mink(x(coneidx), 5)
% if min(x(coneidx)) < 1e-03
%     [~, bid] = mink(x(coneidx), 5);
% else
%     bid = find(x(coneidx) > 1e-05);
% end % End if

bid = find(x(coneidx) > 0);

subcone = m + 1: m + length(bid);
Aid = [(1:m)'; m + bid];

% if length(subcone) ~= ncone
%     keyboard;
% end % End if

if method == "direct"
    d = zeros(n, 1);
    d(coneidx) = x(coneidx).^-2;
    Hess = (- (g * g') / f + ATA) + diag(d * f / rho);
    H11 = Hess(1:m, 1:m); H12 = Hess(1:m, m+1:end);
    H22 = Hess(m + 1:end, m+1:end);
    PH12 = H12 - (H12 * e) * e' / nproj;
    HeeT = (H22 * e) * e' / nproj;
    PH22P = H22 - HeeT - HeeT' + sum(sum(H22)) / nproj^2;
    Hproj = [H11,   PH12;
             PH12', PH22P'];
    [V, evals] = eig(Hproj, 'vector');
    [~, id] = min(evals);
    v = real(V(:, id));
    vbef = [];
elseif method == "sdirect"
    dx = ones(n, 1);
    dx(coneidx) = x(coneidx);
    Xg = dx .* g;
    XATAX = diag(dx) * ATA * diag(dx);
    Hess = rho * (- (Xg * Xg') / f + XATAX);
    H11 = Hess(1:m, 1:m); H12 = Hess(1:m, m+1:end);
    H22 = Hess(m + 1:end, m+1:end) + eye(ncone) * f;
    xe = x(coneidx); nrmxe = norm(xe);
    PH12 = H12 - (H12 * xe) * xe' / nrmxe^2;
    HeeT = (H22 * xe) * xe' / nrmxe^2;
    PH22P = H22 - HeeT - HeeT' + (xe' * H22 * xe) * (xe * xe') / nrmxe^4;
    XHXproj = [H11,   PH12;
               PH12', PH22P'];
    XHXproj = (XHXproj + XHXproj') / 2;
    [V, evals] = eig(XHXproj, 'vector');
    [~, id] = min(evals);
    v = real(V(:, id));
    v = v .* dx;
    vbef = [];
elseif method == "lanczos"
    xsub = x(Aid);
    [v, lam, delta] = potlanczos(xsub, subcone, rho, g(Aid), f, ATA,...
                                 AT(Aid, :), A(:, Aid), false);
    vbef = [];
elseif method == "scaled"
    xsub = x(Aid); 
    [v, lam, delta] = potlanczos(xsub, subcone, rho, g(Aid), f, ATA,...
                                 AT(Aid, :), A(:, Aid), true, vstart);
    vbef = v;
    vbef = vbef / norm(vbef);
    v(subcone) = v(subcone) .* xsub(subcone);
else
    error("Not implemented");
end % End if

vcone = zeros(ncone, 1);
vcone(bid) = v(m + 1:end);
v = [v(1:m); vcone];
v = v / norm(v);

if false
    d = zeros(n, 1);
    d(coneidx) = x(coneidx).^-2;
    Hess = rho * (- (g * g') / f^2 + ATA) + diag(d);
    fprintf("Effect of curvature %f \n", v' * Hess * v);
end % End if 

end % End function