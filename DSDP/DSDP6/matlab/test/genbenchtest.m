clear;
clc;

fname = 'mytest.dat-s';
ns = [3, 4, 5];
K.s = ns;
K.l = 5;

nblock = length(ns);
msps = 5;
mds = 5;
mr1 = 5;

m = msps + mds + mr1;

A = sparse(m, K.l + sum(ns.^2));

for i = 1:msps
    
    A(i, :) = [rand(K.l, 1);
            vec(getrandmat(ns(1))); 
            vec(getrandmat(ns(2)));
            vec(getrandmat(ns(3)))]';  
end % End for

C = [rand(K.l, 1);
     vec(getrandmat(ns(1))); 
     vec(getrandmat(ns(2)));
     vec(getrandmat(ns(3)))]';
 
b = randn(m, 1);

A = full(A);
C = full(C);

writesdpa(fname, A, b, C, K);


function [A] = getrandmat(m) 

x = randi(4);

switch x
    case 1
        A = sprandsym(m, 0.2) + speye(m);
    case 2
        A = randn(m, m);
        A = A' * A;
    case 3
        A = randn(m, 1);
        A = (2 * (randn() > 0.0) - 1) * (A * A');
    case 4
        A = diag(rand(m, 1) * 10);
    case 6
        A = sparse(m, m);
end % End switch

end % End function
    
    
    
    
    