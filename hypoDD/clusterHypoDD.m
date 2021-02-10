function [C,mct,mcc,jE] = clusterHypoDD(mct,mcc,jE,Nobs)
% function [C,mct,mcc,jE] = clusterHypoDD(mct,mcc,Nobs,jE)
% 
% 2020-01-13
% This function divides events into continuously linked clusters for 
% stable double-difference hypocenter inversion
%
%  INPUTS
%
%     mct = [EV1,EV2,STA,PHA,DT,(CC)] pick-based differential-times
%     mcc = [EV1,EV2,STA,PHA,DT,CC] waveform-based differential times
%      jE = Ne x 1 vector of event indices
%    Nobs = minimum number of unique stations/phases to form continuous
%           cluster
%
% OUTPUTS
%   
%      C = Ne x 1 vector of cluster numbers for each event, or 0 if the 
%          event could not be clustered...
%     mct,mcc,jE = (potentially) culled version of inputs

Ne = length(jE);
    
% -- Only need [EV1, EV2, STA, PHA]
% -- In case the event numbers (jE) in mct/mcc aren't 1:Ne, make them so
% -- Ensure events are always reported in same order
A  = [mct(:,1:4); mcc(:,1:4)];
[~,A(:,1:2)] = ismember(A(:,1:2),jE);
A(:,1:2) = sort(A(:,1:2),2);
A = unique(A,'rows');
A = A(:,1:2);

% -- Start counting clusters with negative numbers
C  = (1:Ne)';

% -- Keep iterating until events are fully clustered...
while 1

    % -- Look for cluster-pairs with > Nobs comparisons
    uA = unique(A,'rows');
    [~,juA] = ismember(A,uA,'rows');
    N  = hist(juA,(1:size(uA,1))); 
    jN = find(N >= Nobs);
    
    % -- If no event/cluster pairs have Nobs links, try adding additional
    % -- unclustered events?
    % -- We'd want Nobs station/phase combos at all three events...
    % -- For now, I need at least Nobs comps between a single event pair
    % -- to start a cluster...
    if ~length(jN)
        break
    end
    
    % -- Only change the cluster index of an event once per iteration
    [~,ju] = unique(uA(jN,2));
    jN = jN(sort(ju));
    
    % -- For each pair with >Nobs links... 
    C0 = C;
    for jj = jN 
        C(C==uA(jj,2)) = uA(jj,1);     
    end
    
    % -- Remove comparisons that are now in the same cluster
    [~,jC0] = ismember(A(:,1:2),C0);
    A(:,1:2) = C(jC0);
    A(A(:,1) == A(:,2),:) = [];

    % -- Sort to make sure event/cluster numbers always in same order
    A = sort(A,2);
    
end

% -- Count number of times each cluster index occurs
uC = unique(C);
hC = hist(C,uC);

% -- Remove data with unclustered-events
% -- Remove unclusterd events from 
j0 = uC(find(hC==1));
mct(find(ismember(mct(:,1),jE(j0))+ismember(mct(:,2),jE(j0))),:) = [];
mcc(find(ismember(mcc(:,1),jE(j0))+ismember(mcc(:,2),jE(j0))),:) = [];
jE(j0) = [];
C(j0)  = [];
Ne = length(jE);

% -- Renumber clusters consecutively
[~,C] = ismember(C,unique(C));
Ne = length(jE);
