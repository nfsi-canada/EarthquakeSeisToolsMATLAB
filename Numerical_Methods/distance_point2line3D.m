function D = distance_point2line3D(x,L1,L2)
% function D = distance_point2line3D(x,L1,L2)
%
% 2020-10-15
% This function finds the (minimum) distance between a arbitrary number of
% points and a single line in 3D
% 
%   INPUTS
%
%      x  == Nx3 coordinates of points [X,Y,Z]
%   L1,L2 == length-3 [X Y Z] coordinates of two points on a line
%
%   OUTPUTS
%
%      D == distance from each point to the line, in whatever the input
%           unit was

L1 = L1(:);
L2 = L2(:);
v1 = repmat(L1',size(x,1),1);
v2 = repmat(L2',size(x,1),1);
a = v1-v2;
b = x-v2;

d = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));