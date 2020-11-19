function [hyp,T,picks,stats] = hypoAPPgrid(hyp0,T0,picks,model,E,N,Z,nR,nZ,minZ,maxZ,vpvs)
% function [hyp,T,picks,stats] = hypoAPPgrid(hyp0,T0,picks,model,E,N,Z,nR,nZ,minZ,maxZ,vpvs)
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
%        E = maximum east-west distance (km)
%        N = maximum north-south distance (km)
%        Z = maximum depth difference (km)
%       nR = number of horizontal grid points       
%       nZ = number of depth grid points
%     minZ = minimum allowable depth (km)
%     maxZ = maximum allowable depth (km)
%     vpvs = A constant Vp/Vs applied to whole model 
%       
%
%  OUTPUTS
%
%      hyp = [lat,lon,dep] estimated hypocenter
%        T = estimated origin time (a datetime)
%    picks = [[picks], res (Tpred-Tobs), D epicentral, azimuth, take-off angle]
%    stats = a dictionary with msft,msft0...


% -- Break down the inputs into parts, count number of picks
latE = hyp0(1);
lonE = hyp0(2);
depE = hyp0(3);
latS = picks(:,1);
lonS = picks(:,2);
npk  = size(picks,1);

% -- vfac == a convenient vector to multiply times for S picks by 1.73
vfac = vpvs.^(picks(:,3)-1);

% -- Get relative Observed and Predicted (for hyp0) travel times 
tObs  = picks(:,4);
d     = distance(latE,lonE,latS,lonS,[6371, 0]);
tPred = vfac.*RayTrace(depE,d,model);

% -- Normalize weights?
wght  = npk*picks(:,5)/sum(picks(:,5));
res   = (tObs-tPred)-WeightedMedian(tObs-tPred,wght);
msft  = sum(wght.*abs(res));
res0  = res;
msft0 = msft;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Set up grid

vE = linspace(-E,E,nR)';
vN = linspace(-N,N,nR)';
[DX,DY] = ndgrid(vE,vN);
[latG,lonG] = reckon(latE,lonE,sqrt(DX.^2+DY.^2),90-atan2d(DY,DX),[6371 0]);
latG = latG(:);
lonG = lonG(:);

nG  = length(latG);

% -- Treat depth indepently, because looping over depths is more efficient
% -- Make sure not to extend grid beyon depth bounds
vZ = linspace(max(depE-Z,minZ),min(depE+Z,maxZ),nZ)';

%depG(find((depG<=minZ)+(depG>maxZ))) = [];
%nZ = length(depG);


vfacG = repmat(vfac,1,nG);
tObsG = repmat(tObs,1,nG);
wghtG = repmat(wght,1,nG);

%resG  = zeros(npk,nG);
%resG  = zeros(npk,nG,nZ);
msftG = zeros(nG,nZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Search grid
% -- Loop over grid points (by depth for efficiency)
for iz = 1:nZ
    
    d = distance(repmat(latG',npk,1),repmat(lonG',npk,1),... 
                 repmat(latS,1,nG),repmat(lonS,1,nG),[6371, 0]);    
    
    % -- npk X nG grid... 
    % -- What if I average over a bunch of velocity models??
    tPred = vfacG.*RayTrace(vZ(iz),d,model);
    
    % -- Compute residual and sum for misfit. 
    % -- Should ideally be subtracting a weighted median here (vs. mean),
    % -- but that's somewhat slow, so leaving an option here to use mean...
    resG = (tObsG-tPred)-repmat(WeightedMedian(tObsG-tPred,wghtG),npk,1);
    %resG(:,:,iz) = (tObsG-tPred)-repmat(sum(wghtG0.*(tObsG-tPred)),npk,1);
    
    % -- Sum |residuals| for L1 misfit
    %msftG(:,iz)  =  sum(wghtG.*abs(squeeze(resG(:,:,iz))),1);
    msftG(:,iz)  =  sum(wghtG.*abs(resG),1);
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Assume grid contains all the probability??
% -- Dimensions are [Lat,Lon,Dep]
% -- If picks have different weights what should go in the exponent?
probG = msftG.^(-npk);
probG = reshape(probG,[nR,nR,nZ]);
probGE = sum(sum(probG,3),2);
probGN = sum(sum(probG,3),1)';
probGZ = squeeze(sum(sum(probG,1),2));

% -- Check 95% confidence intervals (it's possible these will just be at 
% -- the bounds of the grid).
hyp = hyp0;
Crng = [1e-5 0.025 0.975 0.9999];

[hE,cE] = WeightedMedian(vE,probGE,Crng);
[hN,cN] = WeightedMedian(vN,probGN,Crng);
[hyp(3),cZ] = WeightedMedian(vZ,probGZ,Crng);

[hyp(1),hyp(2)] = reckon(latE,lonE,sqrt(hE^2+hN^2),90-atan2d(hN,hE),[6371 0]);

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
stats.ENZ   = [cE(4)-cE(1),cN(4)-cN(1),cZ(4)-cZ(1)]/2;



%[latmax,~]  = reckon(latE,lonE,N,0,[6371,0]);
%[latmin,~]  = reckon(latE,lonE,N,180,[6371,0]);
%[~,lonmax]  = reckon(latE,lonE,E,90,[6371,0]);
%[~,lonmin]  = reckon(latE,lonE,E,270,[6371,0]);
%latV = linspace(latmin,latmax,nR);
%lonV = linspace(lonmin,lonmax,nR);
%[latG,lonG] = ndgrid(latV,lonV);


%[hyp(1),cLAT] = WeightedMedian(latV,probGlat,Crng);
%[hyp(2),cLON] = WeightedMedian(lonV,probGlon,Crng);
%[hyp(3),cDEP] = WeightedMedian(depG,probGdep,Crng);

%[DX,DY,depG] = ndgrid(linspace(-E,E,nR),linspace(-N,N,nR),depG);
%[latG,lonG] = reckon(latE,lonE,sqrt(DX.^2+DY.^2),90-atan2d(DY,DX),[6371 0]);
%[hyp(1),bLAT] = WeightedMedian(latG(:),probG(:),Crng);
%[hyp(2),bLON] = WeightedMedian(lonG(:),probG(:),Crng);
%[hyp(3),bDEP] = WeightedMedian(depG(:),probG(:),Crng);

