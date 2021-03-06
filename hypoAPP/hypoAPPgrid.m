function [hyp,T,picks,stats] = hypoAPPgrid(hyp0,T0,picks,model,params)
% function [hyp,T,picks,stats] = hypoAPPgrid(hyp0,T0,picks,model,params)
%
% 2020-03-14
% This function estimates an earthquake hypocenter by grid-seaching around
% an initial hypocenter estimate. For each grid point, travel times are 
% computed with my code RayTrace (based on "ttimes"). 
%
% This does (L1) minimization of: ||(tObs-tPred) - median(tObs-Pred)||_1
%
%
%  2020-07-09
%  instead of finding the very lowest misfit. It finds
%  the weighted median solution, where the "weights" are actual  
%  probabilities
%
%   INPUTS 
%
%     hyp0 = [lat0,lon0,dep0] initial hypocenter
%       T0 = intial origin time (a datetime)
%    picks = [latS,lonS,1/2 for P/S,T (relative to T0 in s), Weight]
%             it may also have columns for [res,dist_epi,azi,toa]
%    model = model name 'GSC','ALASKA','ma2011'
%            OR custom gradient model [Depth (km), Velocity (km s-1)]
%   params = a structure with fields:
%          rENZ = initial maximum distance in [E,N,Z] directions (km)
%            NR = number of horizontal grid points       
%            NZ = number of depth grid points
%            minZ = minimum allowable depth (km)
%            maxZ = maximum allowable depth (km)
%            vpvs = A constant Vp/Vs applied to whole model 
%       
%
%  OUTPUTS
%
%      hyp = [lat,lon,dep] estimated hypocenter
%        T = estimated origin time (a datetime)
%    picks = [[picks], res (Tpred-Tobs), D epicentral, azimuth, take-off angle]
%    stats = a dictionary with fields:
%            res0,res   = initial, final residuals for each pick (s)
%            msft0,msft = sum of absolute values for res0,res
%        bndE,bndN,bndZ = distances [-,+] from final hypocenter that
%                         encompass the probability given in params.CI
%                  rENZ = suggested max. distances [contains 99.99% of probability


% -- Break down the inputs into parts, count number of picks
latS = picks(:,1);
lonS = picks(:,2);
Np   = size(picks,1);

% -- vfac == a convenient vector to multiply times for S picks by 1.73
vfac = params.vpvs.^(picks(:,3)-1);

% -- Get relative Observed and Predicted (for hyp0) travel times 
tObs  = picks(:,4);
d     = distance(hyp0(1),hyp0(2),latS,lonS,[6371 0]);
tPred = vfac.*RayTrace(hyp0(3),d,model);

% -- Normalize weights?
wght  = Np*picks(:,5)/sum(picks(:,5));
res   = (tObs-tPred)-WeightedMedian(tObs-tPred,wght);
msft  = sum(wght.*abs(res));
res0  = res;
msft0 = msft;

% -- Effective degrees of freedom
dof = sum(wght)^2/sum(wght.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%          Set up grid          %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vE = linspace(-params.rENZ(1),params.rENZ(1),params.NH)';
vN = linspace(-params.rENZ(2),params.rENZ(2),params.NH)';
[DX,DY] = ndgrid(vE,vN);
[latG,lonG] = reckon(hyp0(1),hyp0(2),sqrt(DX.^2+DY.^2),90-atan2d(DY,DX),[6371 0]);
latG = latG(:);
lonG = lonG(:);
NG   = length(latG);

% -- Compute epicentral distance from each grid point to each station
dG = distance(repmat(latG',Np,1),repmat(lonG',Np,1),... 
              repmat(latS,1,NG), repmat(lonS,1,NG),[6371, 0]);   

% -- Treat depth indepently, because looping over depths is more efficient
% -- Make sure not to extend grid beyon depth bounds
vZ = linspace(max(hyp0(3)-params.rENZ(3),params.minZ), ... 
              min(hyp0(3)+params.rENZ(3),params.maxZ),params.NZ)';

vfacG = repmat(vfac,1,NG);
tObsG = repmat(tObs,1,NG);
wghtG = repmat(wght,1,NG);
msftG = zeros(NG,params.NZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%         Search the grid         %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Loop over grid points (by depth for efficiency)
for iz = 1:params.NZ
    
    % -- NP X NG grid... 
    tPred = vfacG.*RayTrace(vZ(iz),dG,model);
    
    % -- Compute residual and sum for misfit. 
    % -- Should ideally be subtracting a weighted median here (vs. mean),
    % -- but that's somewhat slow, so leaving an option here to use mean...
    resG = (tObsG-tPred)-repmat(WeightedMedian(tObsG-tPred,wghtG),Np,1);
    %resG = (tObsG-tPred)-repmat(sum(wghtG0.*(tObsG-tPred)),Np,1);
    
    % -- Sum |residuals| for L1 misfit
    msftG(:,iz)  =  sum(wghtG.*abs(resG),1);
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Compute hypocenter as probability centroid       %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Assume grid contains all the probability??
% -- Dimensions are [Lat,Lon,Dep]
% -- If picks have different weights what should go in the exponent?
probG  = msftG.^(-dof);
probG  = reshape(probG,[params.NH,params.NH,params.NZ]);
probGE = sum(sum(probG,3),2);
probGN = sum(sum(probG,3),1)';
probGZ = squeeze(sum(sum(probG,1),2));

% -- Check e.g. 95% confidence intervals 
hyp = hyp0;
CB = [(1-params.CI) (1+params.CI)]/2;
Crng = [5e-5 CB 0.99995];

[hE,cE] = WeightedMedian(vE,probGE,Crng);
[hN,cN] = WeightedMedian(vN,probGN,Crng);
[hyp(3),cZ] = WeightedMedian(vZ,probGZ,Crng);

[hyp(1),hyp(2)] = reckon(hyp0(1),hyp0(2),sqrt(hE^2+hN^2), ... 
                         90-atan2d(hN,hE),[6371 0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     Update the origin and pick times        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Predicted travel times for final hypocenter
[d,azi]     = distance(hyp(1),hyp(2),latS,lonS,[6371, 0]);
[tPred,toa] = RayTrace(hyp(3),d,model);
tPred       = vfac.*tPred;

% -- Compute origin time and adjust picks accordingly
dt0        = WeightedMedian(tObs-tPred,wght);
tObs       = tObs-dt0;
T          = T0+seconds(dt0);

% -- Recompute misfit?
res   = tObs-tPred;
msft  = sum(wght.*abs(res));

picks(:,4) = tObs;
picks      = [picks(:,1:5),res,d,azi,toa];

% -- Build stats dictionary...what else?
stats.msft0 = msft0;
stats.msft  = msft;
stats.res0  = res0;
stats.res   = res; % -- This is also a column in picks
stats.bndE  = cE(2:3)-hE;
stats.bndN  = cN(2:3)-hN;
stats.bndZ  = cZ(2:3)-hyp(3);
stats.rENZ  = [max(abs(cE)),max(abs(cN)),max(abs(cZ-hyp(3)))];
stats.dof   = dof;



