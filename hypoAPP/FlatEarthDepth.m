function ZF = FlatEarthDepth(Z0)
% function ZF = FlatEarthDepth(Z0)
%
% 2020-12-04
% This function converts depths into the earth into "flat earth model"
% depths for a earth radius of 6371 km
%
%   INPUTS
%
%     Z0 == natural depths (km)
%
%  OUTPUTS
%
%     ZF == "flat earth" depths (km)

a  = 6371;
ZF = -a*log((a-Z0)/a);