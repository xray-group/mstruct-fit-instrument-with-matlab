function [dy] = pseudoVoigtAsymDeriv(a,x)
%---------------------------------------------------------------
% PSEUDOVOIGTAsymDERIV    Calculate partial derivatives of
%                         the asymmetric pseudo-Voigt function.
%
% [y] = pseudoVoigtAsymDeriv(a,x)
%
% a - parameters of the asymmetric pseudo-Voigt function
% x - a (row) vector of points where partial derivatives of
%     the asym. pseudo-Voigt function should be calculated.
% 
% y - calculated partial derivatives of the asym. pseudo-Voigt
%     function (5 rows matrix)
%
% example: a derivative of a Gauss-like peak with respect
%          to the position parameter
%
%              x = [-10:0.1:10];
%              a = [1.0, 2.0, 1.5, 0.0, 1.0];
%              dy = pseudoVoigtAsymDeriv(a,x);
%              plot(x,dy(2,:))
%
% For more detailed information see a function 'pseudoVoigtAsym'. 
% 
% In oposite of the function 'pseudoVoigtAsym' this function
% supports only one set of parameters.
%
% X-ray diffraction: see the 'pseudoVoigAsym' function if you
% want to use this function for multiple spectral lines.
%
% See also pseudoVoigtAsym, pseudVoigtAsymFit 
%
% version 1.0, 31.1.2007, (c) Zdenek Matej
%---------------------------------------------------------------
global WAVELENGTHS
% check parameters
if ~exist('WAVELENGTHS','var') | isempty(WAVELENGTHS)
    WAVELENGTHS=[1.0 0.0];
end
if size(x,1)>1  x = reshape(x,1,[]); end
if size(a,2)~=5 a = reshape(a,[],5); end
% calc. partial derivatives (sum over all spectral lines)
dy = zeros(5,length(x));
for k=1:size(WAVELENGTHS,1)
    a2 = a(2)+tan(a(2)*pi/360)*WAVELENGTHS(k,2)*360/pi;
    Lind = x < a(2);
    s1 = (1. +    a(5))*a(3)/2;
    s2 = (1. + 1./a(5))*a(3)/2;
    ss  = (s1 + s2)/2;
    % a1
    c = (1-a(4))*sqrt(log(2)/pi/ss^2);
    dy(1, Lind) = dy(1, Lind) + WAVELENGTHS(k,1)*c*exp(-log(2)/s1^2*(x( Lind)-a2).^2);
    dy(1,~Lind) = dy(1,~Lind) + WAVELENGTHS(k,1)*c*exp(-log(2)/s2^2*(x(~Lind)-a2).^2);
    c = a(4)/pi/ss;
    dy(1, Lind) = dy(1, Lind) + WAVELENGTHS(k,1)*c./(1+(x( Lind)-a2).^2/s1^2);
    dy(1,~Lind) = dy(1,~Lind) + WAVELENGTHS(k,1)*c./(1+(x(~Lind)-a2).^2/s2^2);
    % a2
    c = (1-a(4))*sqrt(log(2)/pi/ss^2);
    c = c*2*log(2)*(x-a2);
    dy(2, Lind) = dy(2, Lind) + WAVELENGTHS(k,1)*c( Lind)/s1^2.*exp(-log(2)/s1^2*(x( Lind)-a2).^2);
    dy(2,~Lind) = dy(2,~Lind) + WAVELENGTHS(k,1)*c(~Lind)/s2^2.*exp(-log(2)/s2^2*(x(~Lind)-a2).^2);
    c = a(4)/pi/ss;
    c = c*2*(x-a2);
    dy(2, Lind) = dy(2, Lind) + WAVELENGTHS(k,1)*c( Lind)/s1^2./(1+(x( Lind)-a2).^2/s1^2).^2;
    dy(2,~Lind) = dy(2,~Lind) + WAVELENGTHS(k,1)*c(~Lind)/s2^2./(1+(x(~Lind)-a2).^2/s2^2).^2;
    % a3
    cc = (1-a(4))*sqrt(log(2)/pi/ss^2);
    c( Lind) = cc*(2*log(2)*(x( Lind)-a2).^2/s1^2-1)/a(3);
    c(~Lind) = cc*(2*log(2)*(x(~Lind)-a2).^2/s2^2-1)/a(3);
    dy(3, Lind) = dy(3, Lind) + WAVELENGTHS(k,1)*c( Lind).*exp(-log(2)/s1^2*(x( Lind)-a2).^2);
    dy(3,~Lind) = dy(3,~Lind) + WAVELENGTHS(k,1)*c(~Lind).*exp(-log(2)/s2^2*(x(~Lind)-a2).^2);
    cc = a(4)/pi/ss;
    c( Lind) = cc*(2*(x( Lind)-a2).^2./(s1^2+(x( Lind)-a2).^2)-1)/a(3);
    c(~Lind) = cc*(2*(x(~Lind)-a2).^2./(s2^2+(x(~Lind)-a2).^2)-1)/a(3);
    dy(3, Lind) = dy(3, Lind) + WAVELENGTHS(k,1)*c( Lind)./(1+(x( Lind)-a2).^2/s1^2);
    dy(3,~Lind) = dy(3,~Lind) + WAVELENGTHS(k,1)*c(~Lind)./(1+(x(~Lind)-a2).^2/s2^2);
    % a4
    c = -sqrt(log(2)/pi/ss^2);
    dy(4, Lind) = dy(4, Lind) + WAVELENGTHS(k,1)*c*exp(-log(2)/s1^2*(x( Lind)-a2).^2);
    dy(4,~Lind) = dy(4,~Lind) + WAVELENGTHS(k,1)*c*exp(-log(2)/s2^2*(x(~Lind)-a2).^2);
    c = 1/pi/ss;
    dy(4, Lind) = dy(4, Lind) + WAVELENGTHS(k,1)*c./(1+(x( Lind)-a2).^2/s1^2);
    dy(4,~Lind) = dy(4,~Lind) + WAVELENGTHS(k,1)*c./(1+(x(~Lind)-a2).^2/s2^2);
    % a5
    dca  = (1-1/a(5)^2)/(2+a(5)+1/a(5));
    dfa1 =  1/(1+a(5));
    dfa2 = -1/(1+a(5))/a(5);
    cc = (1-a(4))*sqrt(log(2)/pi/ss^2);
    c( Lind) = cc*(2*log(2)*(x( Lind)-a2).^2/s1^2*dfa1-dca);
    c(~Lind) = cc*(2*log(2)*(x(~Lind)-a2).^2/s2^2*dfa2-dca);
    dy(5, Lind) = dy(5, Lind) + WAVELENGTHS(k,1)*c( Lind).*exp(-log(2)/s1^2*(x( Lind)-a2).^2);
    dy(5,~Lind) = dy(5,~Lind) + WAVELENGTHS(k,1)*c(~Lind).*exp(-log(2)/s2^2*(x(~Lind)-a2).^2);
    cc = a(4)/pi/ss;
    c( Lind) = cc*(2*(x( Lind)-a2).^2./(s1^2+(x( Lind)-a2).^2)*dfa1-dca);
    c(~Lind) = cc*(2*(x(~Lind)-a2).^2./(s2^2+(x(~Lind)-a2).^2)*dfa2-dca);
    dy(5, Lind) = dy(5, Lind) + WAVELENGTHS(k,1)*c( Lind)./(1+(x( Lind)-a2).^2/s1^2);
    dy(5,~Lind) = dy(5,~Lind) + WAVELENGTHS(k,1)*c(~Lind)./(1+(x(~Lind)-a2).^2/s2^2);
end
dy(2:5,:) = dy(2:5,:)*a(1)/sum(WAVELENGTHS(:,1));
dy(1,:)   = dy(1,:)/sum(WAVELENGTHS(:,1));
return;