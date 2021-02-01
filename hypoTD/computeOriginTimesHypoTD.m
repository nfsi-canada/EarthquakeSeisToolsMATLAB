function T1 = computeOriginTimesHypoTD(T0,H,S,phadir,jE,Ft,Fg,Fv,vpvs)
% function T1 = computeOriginTimesHypoTD(T0,H,S,phadir,jE,Ft,Fg,Fv,vpvs)
%
% 2021-01-15
% This function computes a best-fit origin time after the TD hypocenter
% inversion in hypoTD.
% 
% INPUTS
%
%     T0 == Ne x 1 vector of datetimes
%      H == [x,y,z] hypocenter coordiantes (Ne x 3)
%      S == [x,y,z] station coordinates (Ne x 3)
% phadir == [EV,STA,PHA,T]
%     jE == indices of events that were relocated (<=Ne x 1)
%     Ft == interpolant of travel-time
%     Fg == interpolant of vertical direction cosine
%     Fv == interpolant of velocity at source-depth
%
% OUTPUTS
%
%  T1 == updated origin times (Ne x 1)


T1 = repmat(NaT,length(T0),1);

phadir = phadir(ismember(phadir(:,1),jE),:);
[uES,~,jj] = unique(phadir(:,1:2),'rows');
vfac  = vpvs.^(phadir(:,3)-1);
tPred = RayTrace1DhypoDD(H(uES(:,1),:),S(uES(:,2),:),Ft,Fg,Fv);
tPred = tPred(jj).*vpvs;

for ie = jE'
    jp  = find(phadir(:,1)==ie);
    t0  = phadir(jp,4);
    t1  = tPred(jp);   
    w   = t0.^-1;
    w   = w/sum(w);
    dt0 = WeightedMedian(t0-t1,w);
    
    T1(ie) = T0(ie)+seconds(dt0);
end
    