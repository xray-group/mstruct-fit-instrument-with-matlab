function [params] = CalcParams(x,y)

if isempty(x) x = 1:length(y); end
x = x(:); y = y(:);
[params.maxval, params.nmax] = max(y);
params.A = trapz(x,y);
params.BETA = params.A/params.maxval;
params.M0 = sum(y);
params.M1 = sum(x.*y)/params.M0;
params.M2 = sum(x.^2.*y)/params.M0;
n = [1:length(x)]';
n = find(y <= params.maxval/2 & n <= params.nmax);
n = n(end);
x1 = x(n) + (params.maxval/2 - y(n))/(y(n+1)-y(n))*(x(n+1)-x(n));
n = [1:length(x)]';
n = find(y <= params.maxval/2 & n >= params.nmax);
n = n(1);
x2 = x(n) + (params.maxval/2 - y(n))/(y(n-1)-y(n))*(x(n-1)-x(n));
params.FWHM = x2-x1;
return;