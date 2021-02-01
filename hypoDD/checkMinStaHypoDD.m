function [mct,mcc,jE] = checkMinStaHypoDD(mct,mcc,jE,minsta,minstaPS)
% function [mct,mcc,jE] = checkMinStaHypoDD(mcc,mcc,jE,minsta,minstaPS)
%
% 2021-01-14
% This function counts the numbers of stations with a P comp and the
% number of stations with an S comp for each event, ensuring the total 
% number of unique stations and station/phase combos is met. Ideally 
% we would also make sure the stations are far enough apart...
%
% INPUTS
%
%     mct = [EV1,EV2,STA,PHA,DT,(CC)] pick-based differential-times
%     mcc = [EV1,EV2,STA,PHA,DT,CC] waveform-based differential times
%      jE = vector of event numbers (from 1 to Ne0) still in inversion
%  minsta = minimum number of unique stations per event
% minstaPS= minimum number of unique phases per event, i.e. P and S
%           at the same station are counted separately
%
% OUTPUTS
%  
%      mct,mcc,jE == culled versions of inputs

Ne = length(jE);

while 1

    NstaP = zeros(Ne,1);
    NstaS = zeros(Ne,1);
    Nsta  = zeros(Ne,1);

    MP = [mct(mct(:,4)==1,:); mcc(mcc(:,4)==1,:)];
    MS = [mct(mct(:,4)==2,:); mcc(mcc(:,4)==2,:)];

    for ie = 1:Ne
        staP = unique(MP(find((MP(:,1)==jE(ie))+(MP(:,2)==jE(ie))),3));
        staS = unique(MS(find((MS(:,1)==jE(ie))+(MS(:,2)==jE(ie))),3));
        NstaP(ie) = length(staP);
        NstaS(ie) = length(staS);
        Nsta(ie)  = length(unique([staP(:); staS(:)]));
    end
    
    jcull = find((Nsta <minsta) + (NstaP+NstaS < minstaPS));
    if isempty(jcull)
        break
    else

        % -- Remove equations with culled events
        mct(find(ismember(mct(:,1),jE(jcull)) +  ... 
                 ismember(mct(:,2),jE(jcull))),:) = [];
        mcc(find(ismember(mcc(:,1),jE(jcull)) + ... 
                 ismember(mcc(:,2),jE(jcull))),:) = [];
        jE(jcull) = [];
        Ne = length(jE);

        if Ne < 2
            break
        end
    end
end
    