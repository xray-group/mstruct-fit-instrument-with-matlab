function [a,b,da,db,chi2,Rwp,S,C,n] = pseudoVoigtFit(x,y,w,a0,b0,mu,Linda,Lindb)
%---------------------------------------------------------------
% PSEUDOVOIGTFIT    Fit data by a sum of pseudo-Voigt
%                   functions and a polynomial background
%                   using a simple least squares algorithm.
%
% [a,b,chi2,Rwp,S,da,db,C,n]
%             = pseudoVoigtFit(x,y,w,a0,b0,mu,Linda,Lindb)
%
% x     - x data points (a row vector)
% y     - y data points (a row vector)
% w     - data weights (a row vector)
% a0    - starting parameters of pseudo-Voigt functions,
%         a (n x 4) matrix, where n is the number
%         of pseudo-Voigt functions. For more information
%         see a 'pseudoVoigt' function.
% b0    - a (row) vector specifing the polynomial background,
%         see POLYVAL
% mu    - centering and scaling parameters for the background
%         polynom, see POLYVAL
% Linda - logical matrix (same size as a0) specifing refined
%         pseudo-Voigt functions parameters
% Lindb - logical vector (same size as b0) specifing refined
%         parameters of the background polynomial
% 
% a     - matrix of refined pseudo-Voigt functions parameters
% b     - vector of refined parameters of the background
%         polynomial
% da    - matrix of esds of refined pseudo-Voigt functions
%         parameters
% db    - vector of esds of refined parameters
%         of the background polynomial
% chi2  - chi2 value
% Rwp   - Rwp factor (R-weighted pattern)
% S     - S factor (goodness of fit)
% C     - covariance matrix
% n     - number of iteractions
%
% example: fitting of two peaks
%
%            % create data
%            x = [25:0.05:40];
%            a = [400 33 0.2 0.5; ...
%                 500 31 0.4 0.4];
%            b = [400 -350 200];
%            mu(1) = (x(1)+x(end))/2;
%            mu(2) = x(end)-x(1);
%            y = sum(pseudoVoigt(a,x),1)+polyval(b,x,[],mu);
%            y = y+sqrt(y).*randn(1,length(y));
%            plot(x,y,'k*') 
%            % initial guess
%            a0 = [600 33.15 0.3 0.5; ...
%                  600 30.92 0.3 0.5];
%            Linda = [1 1 1 1;
%                     1 1 1 1];
%            b0 = [0 0 100];
%            Lindb = [1 1 1];
%            yc = sum(pseudoVoigt(a0,x),1)+polyval(b0,x,[],mu);
%            hold on
%            plot(x,yc,'g')
%            % fitting
%            w = 1./y;
%            [a,b] = pseudoVoigtFit(x,y,w,a0,b0,mu,Linda,Lindb)
%            yc = sum(pseudoVoigt(a,x),1)+polyval(b,x,[],mu);
%            plot(x,yc,'b')
%              
% For more detailed information see a function 'pseudoVoigt'. 
%
% X-ray diffraction: see the 'pseudoVoigt' function if you
% want to use this function for multiple spectral lines.
%
% References:
%  [1] W.H.Press, S.A.Teukolsky, W.T.Vetterling, B.P.Flannery,
%       NUMERICAL RECIPES IN C
%
% See also pseudoVoigt, pseudVoigtDeriv 
%
% version 1.0, 14.3.2005, (c) Zdenek Matej
%---------------------------------------------------------------

% check data
if size(x,1)>1 x = reshape(x,1,[]); end
if size(y,1)>1 y = reshape(y,1,[]); end
if size(w,1)>1 w = reshape(w,1,[]); end

if ~exist('mu','var') | isempty(mu) mu = [0.0 1.0]; end

if size(a0,2)~=4 a0 = reshape(a0,[],4); end
if ~exist('Linda','var') | isempty(Linda) Linda = ones(size(a0)); end
if size(Linda,2)~=4 Linda = reshape(Linda,[],4); end    

if ~exist('Lindb','var') | isempty(Lindb) Lindb = ones(size(b0)); end
b0 = b0(:); Lindb = Lindb(:);

Linda = logical(Linda);
Lindb = logical(Lindb);

% initial guess of parameters
a = a0;
b = b0;

na = length(find(Linda));
nb = length(find(Lindb));

% calc scale
yc = sum(pseudoVoigt(a,x),1) + polyval(b,x,[],mu);
s = (w.*yc)*y'/((w.*yc)*yc');

% calc chi2
yc = s*yc;
chi2 = (y-yc).^2*w';
% iteraction cycle
for n=1:100
    % alpha, beta
    D = s*calcD(x,a,b,mu,Linda,Lindb);
    beta = D.*repmat(w,size(D,1),1);
    alpha = beta*D';
    beta  = beta*(y-yc)';
    
    % solve set of equations
    daa = alpha\beta;
    
    % modify parameters
    if na>0
        da = zeros(size(a'));
        da(Linda') = daa(1:na);
        a = a + da';
    end
    if nb>0
        db = zeros(size(b'));
        db(Lindb') = daa(na+1:end);
        b = b + db';
    end
    
    % calc new value chi2
    yc = s*( sum(pseudoVoigt(a,x),1) + polyval(b,x,[],mu) );
    chi20 = chi2;
    chi2 = (y-yc).^2*w';
    
    % terminating condition
    if (chi2<1.e-7 | (chi20-chi2)/chi20<1.e-4) break; end
end

% calc statistic factors and esd
D = s*calcD(x,a,b,mu,Linda,Lindb);

beta = D.*repmat(w,size(D,1),1);
alpha = beta*D';

% covariance matrix
C = inv(alpha);

% esds
daa = sqrt(C(logical(eye(size(alpha)))));
if na>0
    da = zeros(size(a'));
    da(Linda') = daa(1:na);
    da = da';
else
    da = zeros(size(a0));
end
if nb>0
    db = zeros(size(b'));
    db(Lindb') = daa(na+1:end);
else
    db = zeros(size(b0));;
end
if size(b,1)>1 b = reshape(b,1,[]); end

Rwp = sqrt(chi2/(y.^2*w'));
S = sqrt(chi2/(length(x)-na-nb));

% rescale parameters
a(:,1) = s*a(:,1); da(:,1) = s*da(:,1);
b = s*b; db = s*db;

return;

%------------------------------------------------------------
% auxiliary functins
function [D] = calcD(x,a,b,mu,Linda,Lindb)
D = [];
for k=1:size(a,1)
    if sum(Linda(k,:))>0
        d = pseudoVoigtDeriv(a(k,:),x);
        D = [ D; d(Linda(k,:),:) ];
    end
end
for k=1:size(b,1)
    if Lindb(k)
        D = [ D; ((x-mu(1))/mu(2)).^(length(b)-k) ];
    end
end
return;
