function [a,b,ux,uy,alpha,p,chiopt,Cab,Calphap] = robust_wtls(x,y,ux,uy,E,mxiter)
% function [a,b,alpha,p,chiopt,Cab,Calphap] = robust_wtls(x,y,ux,uy,E,mxiter)
%
% INPUTS:
%
%    x,y   = data you are fitting
%    ux,uy = uncertainties of each point 
%    E     = number of standard deviations from mean that is accepted
%    mxiter= maximum number of iterations
%
% robust_wtls finds a total-least squares best-fit line resistant to outliers
% by iteratively increasing the errors of outliers. 
%
k = 0;
N = length(x);
while k < mxiter

    [a,b,alpha,p,chiopt,Cab,Calphap] = wtls_line(x,y,ux,uy);

    % -- Easy formula for distance from line, take stats
    d  = abs(a*x-y+b)./sqrt(a^2+1);

    w  = 1./(ux.*uy);
    md = sum(d.*w)/sum(w);
    dd = N*(abs(d-md).*w)/sum(w);

    ndx = dd/median(dd);

    %wmd = mean(d);
    %sd  = std(d); % of course it will never 'converge' if I don't use some weighted measure
    %nrmd = abs(d-wmd)/sd;

    iout = find(ndx > E);

    if isempty(iout)
        break
    else
        %ux(iout) = ux(iout).*ndx(iout);
        %uy(iout) = uy(iout).*ndx(iout);
        ux(iout) = ux(iout).*(ndx(iout)/E).^2;
        uy(iout) = uy(iout).*(ndx(iout)/E).^2;
        k = k+1;
    end
end
if k >= mxiter
    disp('robust_wtls did not converge!')
end
