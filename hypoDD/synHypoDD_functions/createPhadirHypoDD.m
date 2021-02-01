function [phadir0,phadir1,T1] = createPhadirHypoDD(hyp0,hyp1,T0,sta,model,vpvs,NF0,NF1)
% function [phadir,phadir1,T1] = createPhadirHypoDD(hyp0,hyp1,T0,sta,model,vpvs,NF0,NF1)
% 
% 2021-01-11
% This function takes in hypocenter and station coordiantes and creates 
% a directory of phases...finding travel times in the 3D velocity model.
%
%  INPUTS
%
%    hyp0 = [lat,lon,dep] true hypocenters (Ne x 3)
%    hyp1 = [lat,lon,dep] perturbed/initial hypocenter estimates (Ne x 3)
%      T0 = true origin times (datetimes)  (Ne x 1)
%     sta = [lat,lon,dep]
%   model = a strcuture with fields "mapcmd" and "V"
%    vpvs = (constant) Vp/Vs ratio
%     NF0 = portion of phases to cull e.g. [0.2,0.5] for 20% of P phases
%           and 50% of S phases
%     NF1 = std. dev. in seconds of gaussian noise ot add to picks
%
%  OUTPUTS
%
%     phadir0 = [EV,STA,PHA,T,W] true phase directory
%     phadir1 = [EV,STA,PHA,T,W] initial/perturbed phase directory


col = @(x) x(:);

Ne = size(hyp0,1);
Ns = size(sta,1);

jE = col(repmat(1:Ne,2*Ns,1));
jS = repmat(col(repmat(1:Ns,2,1)),Ne,1);
jP = repmat([1;2],Ne*Ns,1);

phadir0 = [jE,jS,jP,zeros(Ne*Ns*2,1),ones(Ne*Ns*2,1)];
[phadir0,jkeep] = cullPhadirHypoDD(phadir0,NF0);
Np = size(phadir0,1);


[H0x,H0y] = m_ll2xy(hyp0(:,2),hyp0(:,1));
[H1x,H1y] = m_ll2xy(hyp1(:,2),hyp1(:,1));
[Sx,Sy] = m_ll2xy(sta(:,2),sta(:,1));
H0 = [H0x/1000,H0y/1000,hyp0(:,3)];
H1 = [H1x/1000,H1y/1000,hyp1(:,3)];
S  = [Sx/1000,Sy/1000,sta(:,3)];

% -- Solve for travel-times
[upha,~,jj] = unique(phadir0(:,1:2),'rows');
[Ft,Fg,Fv]  = prep1DgridHypoDD(model,10000,500);
t0 = RayTrace1DhypoDD(H0(upha(:,1),:),S(upha(:,2),:),Ft,Fg,Fv);
t1 = RayTrace1DhypoDD(H1(upha(:,1),:),S(upha(:,2),:),Ft,Fg,Fv);

% -- Multiply S-wave times by Vp/Vs             
phadir0(:,4) = t0(jj).*vpvs.^(phadir0(:,3)-1);

% -- Perturb travel times
phadir1 = phadir0;
phadir1(:,4) = phadir1(:,4) + NF1*randn(Np,1);

% -- Perturb origin times according to travel times
T1  = T0;
tp1 = t1(jj).*vpvs.^(phadir1(:,3)-1);


for ie = 1:Ne
    jj = find(phadir1(:,1)==ie);
    dT = mean(phadir1(jj,4)-tp1(jj));
    T1(ie) = T1(ie) + seconds(dT);
    phadir1(jj,4) = phadir1(jj,4)-dT;
end


%dt0    = WeightedMedian(tObs-tPred,wght);
%tObs   = tObs-dt0;
%T      = T+seconds(dt0);


