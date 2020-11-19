function [RT,aziB] = EN2RT(EN,azi,dA,winP)
% function [RT,aziB] = EN2RT(EN,azi,dA,winP)
%
% 2020-06-02
% Rotates east and north components ofa seismogram into radial and 
% tranverse;
%
%   INPUTS
%
%      EN = [east, north] components as columns
%     azi = azimuth to the station in degrees from north (NOT back-azimuth)
%      dA = (optional) a number of degrees to search around azimuth for
%           a "best-fit" rotation, i.e. the azimuth that minimizes 
%           P-energy on the transverse component / maximizes it on the
%           radial component. Always uses 1-degree increments.
%    winP = indices of P-wave window (e.g. capturing 0.5 or 1 sec. around 
%           direct P-arrival)
%
%  OUTPUTS
%
%     RT = [radial, transverse] components as columns
%          away from source is positive on radial
%          90 clockwise from that is positive on transverse
%   aziB = "best=fit" azimuth


if nargin < 3
    a   = 90-azi;
    RT = [cosd(a)*EN(:,1)+sind(a)*EN(:,2), sind(a)*EN(:,1)-cosd(a)*EN(:,2)]; 
    aziB = azi;
else
    dA = fix(dA);
    rngA = azi+[-dA:dA];
    PWR = -Inf;
    for ii = 1:2*dA+1
         a   = 90-rngA(ii);
         rt = [cosd(a)*EN(:,1)+sind(a)*EN(:,2), sind(a)*EN(:,1)-cosd(a)*EN(:,2)]; 
         pwr = sum(rt(winP,1).^2) - sum(rt(winP,2).^2);
         if pwr > PWR
             RT = rt;
             PWR = pwr;
             aziB = rngA(ii);
         end
    end    
end
    


