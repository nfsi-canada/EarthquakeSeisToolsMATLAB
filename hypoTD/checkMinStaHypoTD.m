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

    MP = [mcttd(mcttd(:,5)==1,:); mcctd(mcctd(:,5)==1,:)];
    MS = [mcttd(mcttd(:,5)==2,:); mcctd(mcctd(:,5)==2,:)];

    % -- Count number of P,S stations
    for ie = 1:Ne
        staP = unique(MP(find((MP(:,1)==jE(ie))+(MP(:,2)==jE(ie))),3:4));
        staS = unique(MS(find((MS(:,1)==jE(ie))+(MS(:,2)==jE(ie))),3:4));
        NstaP(ie) = length(staP);
        NstaS(ie) = length(staS);
        Nsta(ie)  = length(unique([staP(:); staS(:)]));     
    end
    
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
    