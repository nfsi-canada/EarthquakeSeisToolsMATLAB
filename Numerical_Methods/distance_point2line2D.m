function D = distance_point2line2D(xp,yp,xl,yl)
% function D = distance_point2line2D(xp,yp,xl,yl)
%
% 2020-10-15
% This function finds the (minimum) distance between a arbitrary number of
% points and a single line in 2D
%
%   INPUTS
%
%      xp,yp == coordinates of points, should be equal length
%      xl,yl == 2x1 vectors of points that describe a line in 2D
%
%   OUTPUTS
%
%      D == distance from each point to the line, in whatever the input
%           unit was

dy = diff(yl);
dx = diff(xl);

D = abs(dy*xp-dx*yp+xl(2)*yl(1)-xl(1)*yl(2))./sqrt(dx^2+dy^2);