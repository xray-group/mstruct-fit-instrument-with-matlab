function [] = fit_lab6_cu
% fit_lab6

% set global wavelength
global Lambda WAVELENGTHS
Lambda = 1.541874; % Kalpha
WAVELENGTHS = [1.0000, (1.540598-Lambda)/Lambda;
               0.5000, (1.544426-Lambda)/Lambda];
         
% read data
td = read('lab6bb-pixcel-ns.xy', 'xy');
xd.x = td(:,1).'; xd.data = td(:,2).'; xd.time = 1; % conversion from matrix to old struct xrdml-format

% read LaB6 PDF hkl list
%hkl_list = read('lab6_cu_pdf2_hkl.txt','gnu');
tt = readtable('lab6_cu_pcd_hkl.csv', 'HeaderLines',1);
% conversion from csv-table to old matrix format and merge overlapping reflections
hkl_list = merge_overlap([tt.(3) tt.(2)*10. tt.(4)]); 

% set the fitting segments
segments = [ 19.0  23.0; ...
             28.0  46.0; ...
             47.5  56.0; ...
             61.0  90.0; ...
             93.0 124.0; ...
            127.0 150.32];

aa = [];
% fit all segments
for s = 1:size(segments,1)
    [a,da] = fit_segment(xd,hkl_list,segments(s,:));
    disp('please check fit and press space to continue')
    pause,
    aa = [aa; a da];
end

fwhm = 2*aa(:,3).*1/4.*(2+aa(:,5)+1./aa(:,5));
figure,plot(aa(:,2),fwhm),xlabel('2Theta (deg)'),ylabel('true-fwhm')
figure,plot(aa(:,2),aa(:,4)),xlabel('2Theta (deg)'),ylabel('shape-k')
figure,plot(aa(:,2),aa(:,5)),xlabel('2Theta (deg)'),ylabel('Asym')

% save data
fid = fopen('params_1.txt','wt');
disp('results saved as params_1.txt'),

fprintf(fid,'#      2Theta (deg)         Intensity (cps)       HWHM (deg)            k            Asym\n');
a = [aa(:,2) aa(:,7) aa(:,1)/xd.time aa(:,6)/xd.time aa(:,3) aa(:,8) aa(:,4) aa(:,9) aa(:,5) aa(:,10)] ;
fprintf(fid,'%10.4f +/- %6.4f   %10.2f +/- %4.2f    %6.4f +/- %6.4f  %4.2f +/- %4.2f  %5.3f +/- %5.3f\n', ...
    a');

fclose(fid);

return,

function [a,da] = fit_segment(xd,hkl_list,segment)
    
    % select fitting interval and prepare data
    Lindx = (xd.x >= segment(1)) & (xd.x <= segment(2));
    x = xd.x(Lindx); y = xd.data(Lindx)*xd.time; w = y; w(w<1) = 1; w = 1./w;
    
    % plot data
    hold off
    plot(x,y,'b*'),
    
    % find peaks in the interval (from hkl_list)
    Lind_hkl = (hkl_list(:,1) >= x(1) ) & (hkl_list(:,1) <= x(end) );
    ind_hkl = find(Lind_hkl,1);
    if isempty(ind_hkl), return, end
   
    % profile fitting
    
    % guess peak params.
    a = prepare_for_fit(x,y,hkl_list(Lind_hkl,:));
    
    % background
    mu = [(x(1)+x(end))/2 x(end)-x(1)];
    b = [0 mean([y(1:4) y(end-5:end)])];
    
    % fix/free fitting parametrs
    Linda = ones(size(a,1),5); %Linda(a(:,2)>60,5) = 0;
    Lindb = [1 1];
    
    ys = sum(pseudoVoigtAsym(a,x),1) + polyval(b,x,[],mu);
    hold on, plot(x,ys,'y-')
    
    [a,b,da,db,chi2,Rwp,S,C,niter] = ...
                pseudoVoigtAsymFit(x,y,w,a,b,mu,Linda,Lindb);
            
    if any( a(:,4)<0.0 | a(:,4)>1.0 )
        warning(['Shape parameter of the pseudo-Voigt function ' ...
                 'out of the range.'])
    end
    if any( a(:,5)<0.1 | a(:,5)>10. )
        warning(['Asymmetry parameter of the pseudo-Voigt function ' ...
                 'out of the range.'])
    end
    
    yc = sum(pseudoVoigtAsym(a,x),1) + polyval(b,x,[],mu);
    hold on, plot(x,yc,'r-',x,y-yc,'g-')
    
    da = S * da;
return,

function [a] = prepare_for_fit(x,y,hkl_list)
    
    global WAVELENGTHS
    
    % peak parameters
    a = zeros(size(hkl_list,1),5);
    
    % for all peaks in the segment
    for k = 1:size(hkl_list,1)
        % select the peak +/- 1 deg range
        Lindx = (x >= hkl_list(k,1)-1) & (x <= hkl_list(k,1)+1);
        xx = x(Lindx); yy = y(Lindx);
        % subtract background
        yy = yy - y(1) + (xx-x(1))/(x(end)-x(1))*(y(end)-y(1));
        % calc peak params.
        p = CalcParams(xx,yy);
        kk = (p.FWHM/2/p.BETA - sqrt(log(2)/pi))/(1/pi-sqrt(log(2)/pi));
        kk = max([kk 0.]); kk = min([kk 1.0]);
        % xx(p.nmax) is probably Kalpha1 - we need correction
        dxx = -2*WAVELENGTHS(1,2)*tan(xx(p.nmax)/360*pi) * 180/pi;
        a(k,:) = [p.A xx(p.nmax)+dxx p.FWHM/2 kk 1.0];
    end
    
return,

function [hkl_out] = merge_overlap(hkl_in)
    hkl_out = zeros(0,size(hkl_in,2));
    for irow=1:size(hkl_in,1)
        idx = find(hkl_out(:,1)==hkl_in(irow,1));
        if isempty(idx)
            hkl_out = vertcat(hkl_out, hkl_in(irow,:));
        else
            % add intensities
            hkl_out(idx,3) = hkl_out(idx,3) + hkl_in(irow,3);
        end
    end
return,
