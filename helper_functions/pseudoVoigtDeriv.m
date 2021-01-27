function [dy] = pseudoVoigtDeriv(a,x)
%-------------------------------------------------------------
% PSEUDOVOIGTDERIV    Calculate partial derivatives of
%                     the pseudo-Voigt function.
%
% [y] = pseudoVoigtDeriv(a,x)
%
% a - parameters of the pseudo-Voigt function
% x - a (row) vector of points where partial derivatives of
%     the pseudo-Voigt function should be calculated.
% 
% y - calculated partial derivatives of the pseudo-Voigt
%     function (4 rows matrix)
%
% example: a derivative of a Gauss-like peak with respect
%          to the position parameter
%
%              x = [-10:0.1:10];
%              a = [1.0, 2.0, 1.5, 0.0];
%              dy = pseudoVoigtDeriv(a,x);
%              plot(x,dy(2,:))
%
% For more detailed information see a function 'pseudoVoigt'. 
% 
% In oposite of the function 'pseudoVoigt' this function
% supports only one set of parameters.
%
% X-ray diffraction: see the 'pseudoVoigt' function if you
% want to use this function for multiple spectral lines.
%
% See also pseudoVoigt, pseudVoigtFit 
%
% version 1.0, 14.3.2005, (c) Zdenek Matej
%-------------------------------------------------------------
global WAVELENGTHS
% check parameters
if ~exist('WAVELENGTHS','var') | isempty(WAVELENGTHS)
    WAVELENGTHS=[1.0 0.0];
end
if size(x,1)>1  x = reshape(x,1,[]); end
if size(a,2)~=4 a = reshape(a,[],4); end
% calc. partial derivatives (sum over all spectral lines)
dy = zeros(4,length(x));
for k=1:size(WAVELENGTHS,1)
    a2 = a(2)+tan(x*pi/360)*WAVELENGTHS(k,2)*360/pi;
    % a1
    c = (1-a(4))*sqrt(log(2)/pi/a(3)^2);
    dy(1,:) = dy(1,:) + WAVELENGTHS(k,1)*c*exp(-log(2)/a(3)^2*(x-a2).^2);
    c = a(4)/pi/a(3);
    dy(1,:) = dy(1,:) + WAVELENGTHS(k,1)*c./(1+(x-a2).^2/a(3)^2);
    % a2
    c = (1-a(4))*sqrt(log(2)/pi/a(3)^2);
    c = c*2*log(2)/a(3)^2*(x-a2);
    dy(2,:) = dy(2,:) + WAVELENGTHS(k,1)*c.*exp(-log(2)/a(3)^2*(x-a2).^2);
    c = a(4)/pi*2*(x-a2)*a(3)./((x-a2).^2+a(3)^2).^2;
    dy(2,:) = dy(2,:) + WAVELENGTHS(k,1)*c;
    % a3
    c = (1-a(4))*sqrt(log(2)/pi/a(3)^2);
    c = c*(2*log(2)*(x-a2).^2-a(3)^2)/a(3)^3;
    dy(3,:) = dy(3,:) + WAVELENGTHS(k,1)*c.*exp(-log(2)/a(3)^2*(x-a2).^2);
    c = a(4)/pi*((x-a2).^2-a(3)^2)./((x-a2).^2+a(3)^2).^2;
    dy(3,:) = dy(3,:) + WAVELENGTHS(k,1)*c;
    % a4
    c = -sqrt(log(2)/pi/a(3)^2);
    dy(4,:) = dy(4,:) + WAVELENGTHS(k,1)*c*exp(-log(2)/a(3)^2*(x-a2).^2);
    c = 1/pi/a(3);
    dy(4,:) = dy(4,:) + WAVELENGTHS(k,1)*c./(1+(x-a2).^2/a(3)^2);
end
dy(2:4,:) = dy(2:4,:)*a(1)/sum(WAVELENGTHS(:,1));
dy(1,:)   = dy(1,:)/sum(WAVELENGTHS(:,1));
return;