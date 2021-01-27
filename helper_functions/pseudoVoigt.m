function [y] = pseudoVoigt(a,x)
%------------------------------------------------------------
% PSEUDOVOIGT    Calculate values of the pseudo-Voigt
%                function.
%
% [y] = pseudoVoigt(a,x)
%
% a - parameters of the pseudo-Voigt function
% x - a (row) vector of points where the pseudo-Voigt
%     should be calculated.
% 
% y - calculated values of the pseudo-Voigt function
%
% example: a Gauss-like peak
%
%              x = [-10:0.1:10];
%              a = [1.0, 2.0, 1.5, 0.0];
%              y = pseudoVoigt(a,x);
%              plot(x,y)
%
% The integral of the pseudo-Voigt function is normalized
% to one and hence the function is defined as:
%
% y = a1 * { (1-a4) * sqrt[log(2)/pi/a3^2]
%                   * exp[-log(2)/a3^2*(x-a2)^2]
%   +        a4  * 1/pi/a3 * 1 / [ 1 + (x-a2)^2/a3^2) ] }
% 
% Parameters have this meaning:
%
% a1 - integrated intensity
% a2 - maximum position
% a3 - HWHM - half width in the half of maximum
% a4 - Cauchy/Gauss component weight (0 - Gauss, 1 - Cauchy)
%
% Function parametes (a) should be optionally stored in
% a vector of the size (1 x 4). If input parameters (a)
% represent a matrix of size (n x 4), n pseudo-Voigt
% functions are returned in n rows of the output vector y.
%
% X-ray diffraction: global WAVELENGTHS
% 
% Function support a global variable WAVELENGTHS that
% represent an array of size (m x 2), where m is a number
% of spectral lines. The WAVELENGTHS matrix has this
% interpretation:
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
% See also pseudoVoigtDeriv, pseudVoigtFit 
%
% version 1.0, 14.3.2005, (c) Zdenek Matej
%------------------------------------------------------------
global WAVELENGTHS
% check parameters
if ~exist('WAVELENGTHS','var') | isempty(WAVELENGTHS)
    WAVELENGTHS=[1.0 0.0];
end
if size(x,1)>1  x = reshape(x,1,[]); end
if size(a,2)~=4 a = reshape(a,[],4); end
% calc. all sets (sum over all spectral lines)
y = zeros(size(a,1),length(x));
for j=1:size(a,1)
    for k=1:size(WAVELENGTHS,1)
        a2 = a(j,2)+tan(a(j,2)*pi/360)*WAVELENGTHS(k,2)*360/pi;
        c = (1-a(j,4))*sqrt(log(2)/pi/a(j,3)^2);
        y(j,:) = y(j,:) + WAVELENGTHS(k,1)*c*exp(-log(2)/a(j,3)^2*(x-a2).^2);
        c = a(j,4)/pi/a(j,3);
        y(j,:) = y(j,:) + WAVELENGTHS(k,1)*c./(1+(x-a2).^2/a(j,3)^2);
    end
    y(j,:) = a(j,1)*y(j,:)/sum(WAVELENGTHS(:,1));
end
return;
