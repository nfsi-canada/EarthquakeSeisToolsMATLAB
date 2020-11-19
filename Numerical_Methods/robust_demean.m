function [mx] = robust_demean(x,w,E,mxiter)
% function [a,b] = robust_demean(x,w,E,mxiter)
%
% INPUTS:
%
%    x     = data you are fitting
%    w     = weights of each point (e.g. cross-corrleation coefficients)
%    E     = number of standard deviations from mean that is accepted
%    mxiter= maximum number of iterations
%
% robust_wtls finds a total-least squares best-fit line resistant to outliers
% by iteratively increasing the errors of outliers. 

k = 0;
N = length(x);
while k < mxiter

    mx  = sum(x.*w)/sum(w);
    dx  = N*(abs(x-mx).*w)/sum(w);
    ndx = dx/median(dx);

    iout = find(ndx > E);

    if isempty(iout)
        break
    else
        w(iout) = w(iout)./(ndx(iout)/E).^2;
        k = k+1;
    end
end
