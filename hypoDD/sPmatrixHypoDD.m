function [GsP,dsP] = sPmatrixHypoTD(msp,H,S,FZsp,w,jE)
% function [GsP,dsP] = sPmatrixHypoTD(msp,H,S,FZsp,w,jE)
%
% 2020-01-15
% This function builds a matrix to enforce sP depth constraits for 
% double-difference hypocenter inversion
% 
% INPUTS
%
%     msp == [EV0,STA,T] 
%       H == [x,y,z] Ne0 x 3 matrix of hypocenters
%       S == [x,y,z] Ne0 x 3 matrix of stations
%            note Ne0 and EV above refer to full (original event list)
%    FZsp == interpolant returning depth given (d,Tsp)
%       w == weight of sP constraints
%      jE == current/retained event indices
%
% OUTPUTS
%
%    GsP == [NsP x 4*Ne] matrix of depth contraints
%    dsP == NsP x 1 vector of (weighted) depths

Ne  = length(jE);
NsP = size(msp,1);
[~,evsSP] = ismember(msp(:,1),jE);

if length(w)==1
    w = repmat(w,NsP,1);
end

% -- Epicentral distance for each sP phase
% -- Check interpolant to get corresponding depth
x  = sqrt(sum((S(msp(:,2),1:2)-H(msp(:,1),1:2)).^2,2));
dz = FZsp(x,msp(:,3))-H(msp(:,1),3);

% -- Form matrix and RHS vector
j1  = (1:NsP)';
j2  = 2*Ne+ evsSP;

% -- Reweight equations to give default equation norms of 1
GsP = sparse(j1,j2,w,NsP,4*Ne);
dsP = w.*dz;

