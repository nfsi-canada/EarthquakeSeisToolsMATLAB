function M = distanceWeightingHypoDD(M,H,dmax)
% function M = distanceWeightingHypoDD(M,H,dmax)
%
% 2021-01-15
% This function culls equations from event pairs greater than a given 
% distance.
%
% INPUTS
%
%      M == mct,mcc,mcttd,or mcctd, matrix of [EVa,EVb,STA...]
%      H == current hypocenters (x,y,z)
%   dmax == maximum distance between events
%
% OUTPUTS
%
%    M == culled version of input M

d  = sqrt(sum((H(M(:,1),:)-H(M(:,2),:)).^2,2));
jj = find(d >= dmax);
M(jj,:) = [];
