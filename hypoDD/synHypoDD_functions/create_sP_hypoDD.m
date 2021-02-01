function msp = create_sP_hypoDD(NsP,errSP,hyp0,hyp1,stadir,phadir,model,vpvs)
% function msp = create_sP_hypoDD(NsP,errSP,hyp0,hyp1,stadir,phadir,model,vpvs)
%
% 2021-01-18
% This creates sP data for synthetic hypoDD/hypoTD problems
%
% INPUTS
%
%      NsP == number of sP phases (will be randomly assigned to events)
%    errSP == standard deviation of error added to sP times (s) 
%     hyp0 == true hypocenters [lat,lon,dep]
%     hyp1 == initial hypocenter estimates [lat,lon,dep]
%   stadir == station coordinates [lat,lon,dep]
%   phadir == true phase directory [EV,STA,PHA,T] (phadir0)
%    model == [DEPTH (km), VEL (km s-1)]
%     vpvs == (constant) Vp/Vs ratio
%
% OUTPUTS
%
%      msp == list of sP phases [EV,STA,TsP]


Ne = size(hyp0,1);

% -- Make sP interpolants
if NsP
    Xmax = 1.5*max(distance(mean(hyp0(:,1)),mean(hyp0(:,2)),stadir(:,1),stadir(:,2),[6371 0]));
    Zmax = 1.5*max([hyp0(:,3); hyp1(:,3)]);
    [FZsp,FTsp] = prepInterpolant_sP(model,vpvs,Xmax,Zmax);
end

% -- Make sP data
msp   = zeros(NsP,3);
msp(:,1) = randi(Ne,[NsP 1]);
for ii = 1:NsP
    phaE = phadir(phadir(:,1)==msp(ii,1),:);
    msp(ii,2) = phaE(randi(size(phaE,1),1),2);
    d_epi     = distance(hyp0(msp(ii,1),1),hyp0(msp(ii,1),2), ... 
                stadir(msp(ii,2),1),stadir(msp(ii,2),2),[6371 0]);
    msp(ii,3) = FTsp(d_epi,hyp0(msp(ii,1),3))+errSP*randn;
end
