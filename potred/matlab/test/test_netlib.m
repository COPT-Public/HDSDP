function [] = test_netlib(fname, maxiter, maxtime, maxmn, minmn, fileID)

data = preprocess(fname);
A = data.A;
b = data.b;
c = data.c;

[m, n] = size(A);

% linesearch = false;
% neweigs = false;

if max(m, n) > maxmn || min(m, n) < minmn
    % Prob pObj dObj pInf dInf rGap Time Status
    fprintf(fileID, "| %30s | %+3.1e | %+3.1e | %3.1e | %3.1e | %3.1e | %5.1f | %s \n",...
        fname, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "Ignored");
    return;
end % End if

params.maxIter = maxiter;
params.maxRuizIter = 0;
params.maxTime = maxtime;
params.coefScal = 0;
params.curvature = 0;
params.curvInterval = 0;
params.RScalFreq = 1000000;
params.PI_RestartMax = 1.0;
params.PI_RestartRate = -1;
params.relFeasTol = 1e-06;
params.relOptTol = 1e-06;

% scaler = max([abs(b); abs(c)]);
% if scaler > 1e+04
%     scaler = 1e+04;
% elseif scaler < 1e-04
%     scaler = 1;
% end % End if

scaler = 1;

tic;

[x, y, s] = potlp(A, b, c, params);
t = toc;

x = x * scaler;
s = s * scaler;
y = y * scaler;

pobj = c' * x;
dobj = b' * y;
pres = norm(A * x - b) / ( 1 + norm(b, 1) );
dres = norm(A' * y + s - c) / ( 1 + norm(c, 1) );
cpl = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));
max1 = max([pres, dres, cpl]);
mm = max1;

if mm < 1e-04
    status = "Optimal";
elseif mm < 1e-03
    status = "Inaccurate";
else
    status = "Failed";
end % End if

fprintf(fileID, "| %30s | %+3.1e | %+3.1e | %3.1e | %3.1e | %3.1e | %5.1f | %s \n",...
        fname, pobj, dobj, pres, dres, cpl, t, status);
end % End function