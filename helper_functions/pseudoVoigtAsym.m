function [y] = pseudoVoigtAsym(a,x)
%--------------------------------------------------------------
% PSEUDOVOIGTASYM    Calculate values of the asymmetric
%                    pseudo-Voigt function.
%
% [y] = pseudoVoigtAsym(a,x)
%
% a - parameters of the asymmetric pseudo-Voigt function
% x - a (row) vector of points where the asymmetric
%     pseudo-Voigt should be calculated.
% 
% y - calculated values of the asym. pseudo-Voigt function
%
% example: a Gauss-like peak
%
%              x = [-10:0.1:10];
%              a = [1.0, 2.0, 1.5, 0.0, 1.0];
%              y = pseudoVoigtAsym(a,x);
%              plot(x,y)
%
% The integral of the asym. pseudo-Voigt function is normalized
% to one and hence the function is defined as:
%
% y = a1 * { (1-a4) * sqrt[log(2)/pi/ss^2]
%                   * exp[-log(2)/s^2*(x-a2)^2]
%   +        a4  * 1/pi/ss * 1 / [ 1 + (x-a2)^2/s^2) ] },
% 
% where ss and s are defined as follows:
%
%                ss = 1/4*(2 + a5 + 1/a5)*a3
%
%         s = 1/2*(1 + a5)*a3          for x <  a2,
%         s = 1/2*(1 + 1/a5)*a3        for x >= a2.
%
% Parameters have this meaning:
%
% a1 - integrated intensity
% a2 - maximum position
% a3 - HWHM - half width in the half of maximum
% a4 - Cauchy/Gauss component weight (0 - Gauss, 1 - Cauchy)
% a5 - asymmetry (1 - symmetric)
%
% Function parametes (a) should be optionally stored in
% a vector of the size (1 x 5). If input parameters (a)
% represent a matrix of size (n x 5), n pseudo-Voigt
% functions are returned in n rows of the output vector y.
%
% X-ray diffraction: global WAVELENGTHS
% 
% Function support a global variable WAVELENGTHS that
% represent an array of size (m x 2), where m is a number
% of spectral lines. (In scripts if for example sigle peaks
% are generated take caution that these variables were not 
% previosly defined by another script.) The WAVELENGTHS matrix
% has this interpretation:
%
% the 1st column - relative line intensity
% the 2st column - delta(Lambda)/Lambda ratio
%                  (Lambda is an averaged wavelength)
%
% example: CuKalpha1 + CuKalpha1 (Lambda = 1.5419)
%
%          WAVELENGTHS =
%              1.0000   -0.0008
%              0.5000    0.0017
% 
% See also pseudoVoigtAsymDeriv, pseudVoigtAsymFit 
%
% version 1.0, 30.1.2007, (c) Zdenek Matej
%--------------------------------------------------------------
global WAVELENGTHS
% check parameters
if ~exist('WAVELENGTHS','var') | isempty(WAVELENGTHS)
    WAVELENGTHS=[1.0 0.0];
end
if size(x,1)>1, x = reshape(x,1,[]); end
if size(a,2)~=5, a = reshape(a,[],5); end
% calc. all sets (sum over all spectral lines)
y = zeros(size(a,1),length(x));
for j=1:size(a,1)
    for k=1:size(WAVELENGTHS,1)
        a2 = a(j,2)+tan(a(j,2)*pi/360)*WAVELENGTHS(k,2)*360/pi;
        Lind = x < a(j,2);
        s1 = (1. +    a(j,5))*a(j,3)/2;
        s2 = (1. + 1./a(j,5))*a(j,3)/2;
        ss  = (s1 + s2)/2;
        c = (1-a(j,4))*sqrt(log(2)/pi/ss^2);
        y(j, Lind) = y(j, Lind) + WAVELENGTHS(k,1)*c*exp(-log(2)/s1^2*(x( Lind)-a2).^2);
        y(j,~Lind) = y(j,~Lind) + WAVELENGTHS(k,1)*c*exp(-log(2)/s2^2*(x(~Lind)-a2).^2);
        c = a(j,4)/pi/ss;
        y(j, Lind) = y(j, Lind) + WAVELENGTHS(k,1)*c./(1+(x( Lind)-a2).^2/s1^2);
        y(j,~Lind) = y(j,~Lind) + WAVELENGTHS(k,1)*c./(1+(x(~Lind)-a2).^2/s2^2);
    end
    y(j,:) = a(j,1)*y(j,:)/sum(WAVELENGTHS(:,1));
end
return;
