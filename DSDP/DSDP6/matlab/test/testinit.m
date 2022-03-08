clear;
clc;

rng(24);

A1 = sprandsym(100, 0.5);
A2 = sprandsym(100, 0.3);
A3 = speye(100);
C = sprandsym(100, 0.5);

Amat = cell(3);
Amat{1} = A1;
Amat{2} = A2;
Amat{3} = A3;
y = zeros(3, 1);

[emin, grad] = eigrad(Amat, C, y);

tic;
for i = 1:100
    y = y + 10 * grad / sqrt(i);
    [emin, grad] = eigrad(Amat, C, y);
    fprintf("%d  %e \n", i, emin);
    if emin > 0.0
        break;
    end % End if
end % End for
toc


opty = 7.355884162036301077147726346084e-02;
optlbd = -1.417483731821011971874213486444e+01;


cvx_begin sdp
variable x(3, 1)
maximize lambda_min(C - x(1) * A1 - x(2) * A2 - x(3) * A3)
cvx_end
