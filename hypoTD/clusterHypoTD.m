function [C,mcttd,mcctd,jE] = clusterHypoTD(mcttd,mcctd,jE,Nobs)
% function [C,mcttd,mcctd,jE] = clusterHypoTD(mcttd,mcctd,Nobs,jE)
% 
% 2020-01-13
% This function divides events into continuously linked clusters for 
% stable double-difference hypocenter inversion
%
%  INPUTS
%
%   mcttd = [EVa,EVv,STA1,STA2,PHA,DT,(CC)] pick-based differential-times
%   mcctd = [EVa,EVb,STA1,STA2,PHA,DT,CC] waveform-based differential times
%      jE = Ne x 1 vector of event indices
%    Nobs = minimum number of unique stations/phases to form continuous
%           cluster
%
% OUTPUTS
%   
%      C = Ne x 1 vector of cluster numbers for each event, or 0 if the 
%          event could not be clustered...
%   mcttd,mcctd,jE = (potentially) culled version of inputs

Ne = length(jE);
    
% -- Only need [EV1, EV2, STA1, STA2, PHA]
% -- In case the event numbers (jE) in mct/mcc aren't 1:Ne, make them so
% -- Ensure events are always reported in same order
A  = [mcttd(:,1:5); mcctd(:,1:5)];
[~,A(:,1:2)] = ismember(A(:,1:2),jE);
A(:,1:2) = sort(A(:,1:2),2);
A = unique(A,'rows');

% -- Start counting clusters with negative numbers
C  = zeros(Ne,1);
Nc = 0;

% -- Keep iterating until events are fully clustered...
while 1

    % -- Look for event-pair with most comparisons
    uEP = unique(A(:,1:2),'rows');
    [~,juEP] = ismember(A(:,1:2),uEP,'rows');
    NEP = hist(juEP,(1:size(uEP,1)));
    
    jN = find(NEP >= Nobs);
    
    % -- If no event/cluster pairs have Nobs links, try adding additional
    % -- unclustered events?
    % -- We'd want Nobs station/phase combos at all three events...
    % -- For now, I need at least Nobs comps between a single event pair
    % -- to start a cluster...
    if ~length(jN)
        break
    end
    
    % -- For each pair with >Nobs links...
    for jj = jN

        evs = uEP(jj,:)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Adding to existing cluster
        if sum(evs < 0) == 1
            jC = -evs(evs<0);
            evN = evs(evs>0);
            C(evN) = jC;
            A(A(:,1:2)==evN) = -jC;
            uEP(uEP==evN) = -jC;
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Joining two clusters    
        elseif sum(evs < 0) == 2
            if evs(1)==evs(2)
                continue
            end
            jC = -evs;
            C0 = min(jC);
            C1 = max(jC);
            C(C==C1) = C0;
            C(C> C1) = C(C > C1)-1;
            A(A(:,1:2)==-C1) = -C0;
            A(A(:,1:2)< -C1) = A(A(:,1:2)< -C1)+1;
            uEP(ismember(uEP,evs)) = -C0;
            
            Nc = Nc-1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Forming a new cluster        
        else
            Nc = Nc + 1;
            C(evs) = Nc;
            A(ismember(A(:,1:2),evs)) = -Nc;
            uEP(ismember(uEP,evs)) = -Nc;
        end
       
        A(find((A(:,1)<0).*(A(:,2)<0)),:) = [];
    end
  
    % -- Sort to make sure event/cluster numbers always in same order
    A(:,1:2) = sort(A(:,1:2),2);
    
end

% -- Remove data with unclustered-events
% -- Remove unclusterd events from 
j0 = find(C==0);
mcttd(find(ismember(mcttd(:,1),jE(j0))+ismember(mcttd(:,2),jE(j0))),:) = [];
mcctd(find(ismember(mcctd(:,1),jE(j0))+ismember(mcctd(:,2),jE(j0))),:) = [];
jE(j0) = [];
C(j0)  = [];
Ne = length(jE);
