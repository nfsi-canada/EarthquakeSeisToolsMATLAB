function [PSV,thtB] = RU2PSV(RU,tht,a,b,d_tht,winP,winS)
% function PSV,thtB] = RU2PSV(RU,tht,a,d_tht,winP,winS)
%
% 2020-04-01
% This function converts radial and vertical (up) components of a 
% seismogram to SV and P components, given a ray parameter and P and S 
% velocities at the surface (a and b). The equation comes from Shearer
% (Introduction to Seismology) 7.34 (on page 200);
%
% THIS CODE ASSUMES THE UPCOMING WAVE IS A P-WAVE
%
%  INPUTS
%
%     RU = [radial component, up component]
%    tht = incidence angle in degrees, 0 = vertically up
%      a = P velocity at the surface
%      b = S velocity at the surface
%  d_tht = (optional) a number of degrees to search around incidence angle
%           for a "best-fit" rotation, i.e. the "tht" that minimizes 
%           P-energy on the P-component / minimizes it on the SV component
%           / minimizes S-energy on the P-component / maximizes it on the
%           S-component. Always uses 1-degree increments.
%   winP = indices of P-wave window (e.g. capturing 0.5 or 1 sec. around 
%          direct P-arrival)
%   winS = indices of S-wave window (can be different length than winP)
%
%  OUTPUTS
%
%    PSV = [P component, SV component]
%   thtB = "best-fit" incidence angle

if nargin < 5 
    p = sind(tht)/a;
    Na = sqrt(a^-2 - p^2);
    Nb = sqrt(b^-2 - p^2);

    % -- X = [X_pr X_pu; X_sr X_su]
    %X = [p*b^2/a,               (1-2*b^2*p^2)/(2*a*Na)
    %    (1-2*b^2*p^2)/(2*b*Nb), -p*b];
    %
    %PSV = (X*RU')';
    
    % -- [P SV] = [R U][X_PR  X_SR
    %                   X_PU  X_SU] 
    X = [pa*b^2/a, (1-2*b^2*pb^2)/(2*b*Nb)         
        (1-2*b^2*pa^2)/(2*a*Na),    -pb*b];
    PSV = RU*X;
    
    thtB = tht;
else
    d_tht = fix(d_tht);
    rngT  = [max(tht-d_tht,1):min(tht+d_tht,89)];
    PWR   = -Inf;
    
    for ii = 1:length(rngT)
        p = sind(rngT(ii))/a;
        Na = sqrt(a^-2 - p^2);
        Nb = sqrt(b^-2 - p^2);

        % -- X = [X_PR X_SR; X_PU X_SU]
        X = [pa*b^2/a, (1-2*b^2*pb^2)/(2*b*Nb)         
            (1-2*b^2*pa^2)/(2*a*Na),    -pb*b];
        psv = RU*X;

        pwr = sum(psv(winP,1).^2) + sum(psv(winS,2).^2) - ...
              sum(psv(winS,1).^2) - sum(psv(winP,2).^2);

        if pwr > PWR
            PSV  = psv;
            PWR  = pwr;
            thtB = rngT(ii);
        end
    end
end
            
%SvP = flipud(SvP)';

