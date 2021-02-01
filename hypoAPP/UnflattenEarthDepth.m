function Z0 = FlatEarthDepth(ZF)
% function Z0 = FlatEarthDepth(ZF)
%
% 2020-12-04
% This function converts depths from "flat earth model" depths to natural
% depths for a earth radius of 6371 km
%
%   INPUTS
%
%     ZF == "flat earth" depths (km)
%
%  OUTPUTS
%
%     Z0 == natural depths (km)

a  = 6371;
Z0 = a*(1-exp(-ZF/a));