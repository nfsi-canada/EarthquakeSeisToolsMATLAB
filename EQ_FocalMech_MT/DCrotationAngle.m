function phi = DCrotationAngle(sdr1,sdr2)
% function phi = DCrotationAngle(sdr1,sdr2)
%
% 2020-07-06
% This function computes the minimum rotation angle between pairs of
% double-couple focal mechanisms. This follows the algorithm of 
% Kagan [2007,GJI]
%
%   INPUTS
%
%      sdr1 = [strike,dip,rake] 
%      sdr2 = [strike,dip,rake] (must be same height as sdr1)
%
%  OUTPUTS
%
%       phi = (minimum) rotation angle, in degrees

N = size(sdr1,1);

% -- Compute Pressure and Tension axes, then Intermediate axes
[~,~,p1,t1] = sdr2PTaxes(sdr1);
[~,~,p2,t2] = sdr2PTaxes(sdr2);
b1 = cross(p1,t1);
b2 = cross(p2,t2);

% -- First assume phi < 90, take absolute value of all dot products 
ppttbb = abs([dot(p1,p2,2),dot(t1,t2,2),dot(b1,b2,2)]);
phi    = acosd(0.5*(sum(ppttbb,2)-1));

% -- Next assume phi > 90, take smallest |dot product| to be negative
[~,jmin]     = min(ppttbb,[],2);
jmin         = sub2ind([N,3],[1:N]',jmin);
ppttbb(jmin) = -ppttbb(jmin);
phi2         = acosd(0.5*(sum(ppttbb,2)-1));

% -- Check when phi > 90
phi(phi >= 90) = phi2(phi >= 90); 


