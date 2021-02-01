function w = distanceWeightingHypoDD(M,H,dmax)
% function w = distanceWeightingHypoDD(M,H,dmax)
%
% 2021-01-15
% This function weights equations based on inter-event distance. Note that
% I follow the cubic weighting sugested by Waldhauser and Ellsworth
% (2000,BSSA) except that no equation is weighted below 0.1.
%
%
% INPUTS
%
%       M == mct,mcc,mcttd,or mcctd, matrix of [EVa,EVb,STA...]
%       H == current hypocenters (x,y,z)
%    dmax == maximum inter-event distance 
%
% OUTPUTS
%
%    w == a weight for each equation in M;

d     = sqrt(sum((H(M(:,1),:)-H(M(:,2),:)).^2,2));
w     = 1-0.9*(d/dmax).^3;
j1    = find(d > dmax);
w(j1) = 0.1;



% I also demean the weights so that the relative weightes of different 
% data types can be retained. If there are interevent distances greater
% they dmax, they get set weights of 0.1
%w  = w/mean(w);




%j1 = find((d > d0).*(d < d1+d0));
%j2 = find(d >= d1+d0);
%w = 1-0.1*(d/d0);
%w(j1) = 0.9-0.8*((d(j1)-d0)/d1);
%w(j2) = 0.1;
%w = w/mean(w);
