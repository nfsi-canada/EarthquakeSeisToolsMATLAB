function [sdr1,sdr2] = PTaxes2sdr(P,T)
% function [sdr1,sdr2] = PTaxes2sdr(P,T)
%
% 2020-07-20
% This function takes in a set of Pressure and Tension axes and returns
% the corresponding set of strike/dip/rake angles.
%
%   INPUTS
%
%       P == Pressure axis. As [x,y,z] or [trend,plunge] in degrees
%       T == Tension axis. As [x,y,z] or [trend,plunge] in degrees
%
%   OUTPUTS
%
%     sdr1,sdr2 == equivalent strike/dip/rake focal mechanisms

if length(P)==2
    P = [cosd(90-P(1))*cosd(P(2)),sind(90-P(1))*cosd(P(2)),-sind(P(2))];
    T = [cosd(90-T(1))*cosd(T(2)),sind(90-T(1))*cosd(T(2)),-sind(T(2))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- First nodal plane
n = (T+P)/sqrt(2);
v = (T-P)/sqrt(2);

% -- Ensure normal is pointing upward 
if n(3) < 0
    n = -n;
    v = -v;  
end

% -- Dip can now be taken without complication
dip = acosd(n(3));

% -- Horizontal angle of normal (dip direction) 
% -- from -180 to 180 (W to W) counter-clockwise
str = -atan2d(n(2),n(1));
str(str<0) = str+360;

rke  = asind(v(3)/sind(dip));
cosR = v(1)*cosd(90-str)+v(2)*sind(90-str);
rke(cosR<0) = sign(rke)*180-rke;


sdr1 = [str,dip,rke];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Second nodal plane
v = (T+P)/sqrt(2);
n = (T-P)/sqrt(2);

% -- Ensure normal is pointing upward 
if n(3) < 0
    n = -n;
    v = -v;  
end

% -- Dip can now be taken without complication
dip = acosd(n(3));

% -- Horizontal angle of normal (dip direction) 
% -- from -180 to 180 (W to W) counter-clockwise
str = -atan2d(n(2),n(1));
str(str<0) = str+360;

rke  = asind(v(3)/sind(dip));
cosR = v(1)*cosd(90-str)+v(2)*sind(90-str);
rke(cosR<0) = sign(rke)*180-rke;

sdr2 = [str,dip,rke];
