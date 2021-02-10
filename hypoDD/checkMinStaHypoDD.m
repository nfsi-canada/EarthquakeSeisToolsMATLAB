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
   
    MP = [mct(mct(:,4)==1,1:3); mcc(mcc(:,4)==1,1:3)];
    MS = [mct(mct(:,4)==2,1:3); mcc(mcc(:,4)==2,1:3)];

    % -- Count number of P,S stations
    AP    = unique([MP(:,[1 3]); MP(:,[2 3])],'rows');
    AS    = unique([MS(:,[1 3]); MS(:,[2 3])],'rows');
    APS   = unique([AP; AS],'rows');
    NstaP = hist(AP(:,1),jE);
    NstaS = hist(AS(:,1),jE);
    Nsta  = hist(APS(:,1),jE);
    
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
    