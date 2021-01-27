function [af,daf] = pseudoVoigtDeconvN(ah,ag,dah)
%ah = [1 0 0.5 1]; dah = [0 0 0.1 0.4];
%ag = [1 0 0.05 1];
%
af = pseudoVoigtDeconvNproc(ah,ag);
if exist('dah','var') & ~isempty(dah)
    dwhm = zeros(1,2); dk = zeros(1,2);
    % hwhm
    if abs(dah(3))>eps
        dd = abs(dah(3)/2);
        a = ah; a(3)=a(3)+dd; 
        a = pseudoVoigtDeconvNproc(a,ag);
        dhwhm(1) = (a(3)-af(3))/dd;
        dhwhm(2) = (a(4)-af(4))/dd;
    end
    % k
    if abs(dah(4))>eps
        dd = abs(dah(4)/2);
        if ah(4)-dd<0.0 & ah(4)+dd>1.0
            warning(['Shape parameter of the pseudo-Voigt function ' ...
                     'out of range.'])
        end
        if ah(4)+dd>1.0 & ah(4)-dd>=0.0 dd=-dd; end
        a = ah; a(4)=a(4)+dd;
        a = pseudoVoigtDeconvNproc(a,ag);
        dk(1) = (a(3)-af(3))/dd;
        dk(2) = (a(4)-af(4))/dd;
    end
    daf = zeros(1,4);
    daf(3) = abs(dhwhm(1))*dah(3)+abs(dk(1))*dah(4);
    daf(4) = abs(dhwhm(2))*dah(3)+abs(dk(2))*dah(4);
end
return;

function [af] = pseudoVoigtDeconvNproc(ah,ag)
epsq = 1e-4; nq = 100;
%
xmax = 0.5/ag(3)*nq;
N = nq/sqrt(epsq)*ah(3)/ag(3); N = ceil(log(N)/log(2)); N = 2^N;
%
x  = linspace(0,xmax,N);
fh = (1-ah(4))*exp(-pi^2*ah(3)^2/log(2)*x.^2) ...
   + ah(4)*exp(-2*pi*ah(3)*abs(x));
fg = (1-ag(4))*exp(-pi^2*ag(3)^2/log(2)*x.^2) ...
   + ag(4)*exp(-2*pi*ag(3)*abs(x));
ff = fh./(fg+epsq);
%semilogy(x,fh,'r',x,fg,'b',x,ff,'g')
ff = [ff fliplr(ff(1:end))];
f  = real(fftshift(ifft(ff,2*N)))*2*xmax;% f(f<0) = eps;
q  = 1/(x(2)-x(1))*linspace(-N+1,N,2*N)/2/N;
p = CalcParams(q,f);
A = trapz(q,f);
hwhm = p.FWHM/2;
beta = A/f(N);
k = (hwhm/beta-sqrt(log(2)/pi))/(1/pi-sqrt(log(2)/pi));
af = [ah(1)/ag(1) 0.0 hwhm k];
fpv = (1-af(4))*sqrt(log(2)/pi/af(3)^2)*exp(-log(2)*q.^2/af(3)^2) ...
    + af(4)/pi/af(3)./(1+q.^2/af(3)^2);
%figure, semilogy(q,f,'*',q,fpv,'r')
return;