function [n,v] = sdr2nv(sdr)
% function [n,v] = sdr2nv(sdr)
%
% 2019-11-27
% Fault normals and slip vectors are computed given strike,dip,rake angles
% This guarantees normal will be pointing up.
% 
%      INPUTS
%
%     sdr = NEx3 matrix of strike,dip,rake angles
%           strike: 0 to 360 in degrees E of N
%           dip:    0 to 90 (right-hand rule from strike)
%           rake:   -180 to +180, where 0 is left-lateral strike-slip
%                                       +90 is reverse
%                                       -90 is normal
%      OUTPUTS
%
%     n = NEx3 matrix of fault normals in ENU convetion
%     v = NEx3 matrix of slip vectors (tracing motion of hanging wall)
%         in ENU convention
%

s = sdr(:,1);
d = sdr(:,2);
r = sdr(:,3);

n = [sind(d).*cosd(-s), sind(d).*sind(-s), cosd(d)];

v = [cosd(r).*sind(s) - cosd(d).*sind(r).*cosd(s), ...
     cosd(r).*cosd(s) + cosd(d).*sind(r).*sind(s), ... 
     sind(d).*sind(r)];

