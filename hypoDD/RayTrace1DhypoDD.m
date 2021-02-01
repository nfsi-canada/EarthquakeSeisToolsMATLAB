function [t,g,vto] = RayTrace1DhypoDD(H,S,Ft,Fg,Fv);
% function [t,g,vto] = RayTrace1DhypoDD(H,S,Ft,Fg,Fv);
%
% 2021-01-14 
% This takes function takes in a hypocenter and station coordinates
% and outputs travel-times, take-off directions, and source-velocities
% using interpolants created by "prep1DgridHypoDD".
%
% INPUTS
%
%     H = [x,y,z] Nray x 3 of hypocenters 
%     S = [x,y,z] Nray x 3 of stations
%    Ft = interpolant of travel-time (d_epi(km),z_eff(km))
%    Fg = interpolant of take-off angle (0=up) (d_epi(km),z_eff(km))
%    Fv = interpolant of source-velocity (z_eff(km))
%
% OUTPUTS
%
%    tpp == predicted travel time (in units matching V)
%      g == take-off directions (Nray x 3 matrix), 
%     v0 == velocity at sources

% -- Compute epicentral distance, effective depth
dx = S(:,1)-H(:,1);
dy = S(:,2)-H(:,2);
d  = sqrt(dx.^2+dy.^2);
z  = H(:,3)-S(:,3);

% -- Interpolate times
t = Ft(d,z);

% -- Interpolate take-off angles
gz = Fg(d,z);
gH = sqrt(1-gz.^2);
g  = [(dx./d).*gH, (dy./d).*gH, gz];

% -- Interpolate source velocities
vto  = Fv(z);