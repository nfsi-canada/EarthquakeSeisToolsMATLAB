function [mcttd,mcctd,jE] = checkMinStaHypoTD(mcttd,mcctd,jE,minsta,minstaPS)
% function [mcttd,mcctd,jE] = checkMinStaHypoTD(mcctd,mcctd,jE,minsta,minstaPS)
%
% 2021-01-15
% This function counts the numbers of stations with a P comp and the
% number of stations with an S comp for each event, ensuring the total 
% number of unique stations and station/phase combos is met. Ideally 
% we would also make sure the stations are far enough apart...
%
% Note that this was modifed from "checkMinStaHypoDD"
%
% INPUTS
%
%   mcttd = [EVa,EVv,STA1,STA2,PHA,DT,(CC)] pick-based differential-times
%   mcctd = [EVa,EVb,STA1,STA2,PHA,DT,CC] waveform-based differential times
%      jE = vector of event numbers (from 1 to Ne0) still in inversion
%  minsta = minimum number of unique stations per event
% minstaPS= minimum number of unique phases per event, i.e. P and S
%           at the same station are counted separately
%
% OUTPUTS
%  
%    mcttd,mcctd,jE == culled versions of inputs

Ne = length(jE);

while 1

    NstaP = zeros(Ne,1);
    NstaS = zeros(Ne,1);
    Nsta  = zeros(Ne,1);

    MP = [mcttd(mcttd(:,5)==1,1:4); mcctd(mcctd(:,5)==1,1:4)];
    MS = [mcttd(mcttd(:,5)==2,1:4); mcctd(mcctd(:,5)==2,1:4)];

    % -- Count number of P,S stations
    AP    = unique([MP(:,[1 3]); MP(:,[1 4]); MP(:,[2 3]); MP(:,[2 4])],'rows');
    AS    = unique([MS(:,[1 3]); MS(:,[1 4]); MS(:,[2 3]); MS(:,[2 4])],'rows');
    APS   = unique([AP; AS],'rows');
    NstaP = hist(AP(:,1),jE);
    NstaS = hist(AS(:,1),jE);
    Nsta  = hist(APS(:,1),jE);

    
    % -- Check for events that don't meet minimum requirements 
    jcull = find((Nsta <minsta) + (NstaP+NstaS < minstaPS));
    if isempty(jcull)
        break
    else
        % -- Remove equations with culled events
        mcttd(find(ismember(mcttd(:,1),jE(jcull)) +  ... 
                   ismember(mcttd(:,2),jE(jcull))),:) = [];
        mcctd(find(ismember(mcctd(:,1),jE(jcull)) + ... 
                   ismember(mcctd(:,2),jE(jcull))),:) = [];
        jE(jcull) = [];
        Ne = length(jE);

        if Ne < 2
            break
        end
    end
    
end
    

%staP = unique(MP(find((MP(:,1)==jE(ie))+(MP(:,2)==jE(ie))),3:4));
%staS = unique(MS(find((MS(:,1)==jE(ie))+(MS(:,2)==jE(ie))),3:4));