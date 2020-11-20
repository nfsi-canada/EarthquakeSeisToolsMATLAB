function [AFM,bFM,NFM] = build_AFM(PFM)
% function [AFM,bFM,NFM] = build_AFM(PFM)
%
% 2020-11-19
% This function makes a linear system out of P-wave first-motion data,
% where the solution vector would contain unit moment tensor elements.
%
%  INPUTS
%
%    PFM = [AZI,TOA,+1/-1,SNR]
%            AZI == azimuth, degrees east of north
%            TOA == take-off angle in degrees
%                     0 == up 
%                    90 == horizontal
%                   180 == down
%          +1/-1 == either +1 for up-first motion or -1 for down-first
%            SNR == signal-to-noise ratio
%
%  OUTPUTS
%
%   AFM = an NFM x 6 matrix of coefficients corresponding to the 
%         moment tensor elements [M11 M22 M33 M12 M13 M23]
%   bFM = NFM x 1 right-hand side vector of +1 or -1
%   NFM = number of P-wave first-motions in PFM

% -- Compute direction cosines of P-waves (as they leave focal sphere)
gx  = cosd(90-PFM(:,1)).*sind(PFM(:,2));
gy  = sind(90-PFM(:,1)).*sind(PFM(:,2));
gz  = cosd(PFM(:,2));

% -- Convert to linear system corresponding to moment tensor elements
% -- [M11 M22 M33 M12 M13 M23]
AFM = [gx.^2, gy.^2, gz.^2, 2*gx.*gy, 2*gx.*gz, 2*gy.*gz];

% -- Right-hand side vector is simply +1/-1
bFM = PFM(:,3);
NFM = size(PFM,1);