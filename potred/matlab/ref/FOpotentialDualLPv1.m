% This is one-iteration implemetation of the first-order potential 
% reduction algorithm for solving
%     max    b^Ty     s.t. A^Ty+s=c, s\ge 0
% The algorithm apply the steepest ascent method to reduce potential 
% function 
%   \rho\ln(z-b'*y) -\sum \ln(s_j), s.t. A'*y+s=c
% where z is a upper bound on b^Ty and (y,s) is interior feasible,
% or
%   \rho\ln(x'*s) -\sum \ln(s_j), s.t. A'*y+s=c, 
% where x is a primal feasible solution and updated in the process.
%
% Set initial data
%
if (iter == 0),
    [m,n]  = size(A);
    ee=ones(n,1);
    s=ee;
    c=ee;
    x=rx;
    b=A*x;
    y=zeros(m,1);
    gap=x'*s
    AAT=inv(A*A');
    rho=n+sqrt(n);
    alpha=0.5;
    toler=1.e-5;
end;
%
while (iter < 500),
%
% compute the gradient of potential function
%
invs = 1./s;
gs = (rho/gap)*x - invs;
%
% compute the gradient projection
%
dy1 = AAT*b;
dy2 = AAT*(A*invs);
x1 = A'*dy1;
x2 = invs-A'*dy2;
%
% compute a possibly improved primal
%
xx =x1+(gap/rho)*x2;
if (xx >= 0)&(s'*xx < gap),
    gap = s'*x1/(1-s'*x2/rho);
    if (x1+(gap/rho)*x2 >= 0),
       x = x1+(gap/rho)*x2;
    else
       theta=max(-x1./max(toler^2,x2));
       x = x1+theta*x2;
       gap=s'*x;
    end;
end;
% 
% Finalize the dual directions
%
dy = (rho/gap)*dy1-dy2;
ds = -A'*dy;
%
% compute the stepsize
%
beta = alpha/norm(ds./s);
%
% scale the projected gradient
%
y = y + beta*dy;
s = s + beta*ds;
iter = iter+1;
gap = x'*s;
end
gap
%

