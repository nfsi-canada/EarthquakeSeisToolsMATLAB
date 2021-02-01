function [hyp,T,picks,stats] = hypoAPP(hyp,T,picks,model,params)
% function [hyp,T,picks,stats] = hypoAPP(hyp,T,picks,model,params)
%
% 2020-03-14
% This function estimates an earthquake hypocenter by grid-seaching around
% an initial hypocenter estimate (calling the function hypoAPPgrid). It
% repeats this process
%
%
%
%
%   INPUTS 
%
%     hyp0 = [lat0,lon0,dep0] initial hypocenter
%       T0 = intial origin time (a datetime)
%    picks = [latS,lonS,1/2 for P/S,T (relative to T0 in s), Weight]
%             it may also have columns for [res,dist_epi,azi,toa]
%    model = model name 'GSC','ALASKA','ma2011'
%            OR custom gradient model [Depth (km), Velocity (km)]
%   params = an (optional) structure with the following fields:
%         rENZ = initial maximum distance in [E,N,Z] directions (km)
%           NR = number of horizontal grid points 
%           NZ = number of depth grid points
%                It makes sense for nR,nZ to be odd numbers so that hyp0
%                is one of the grid points
%         minZ = minimum allowable depth (km)
%         maxZ = maximum allowable depth (km)
%          tol = stop iterating when grid size less than "tol" km
%      picktol = cull picks with residual > picktol * std. dev.
%         vpvs = constant Vp/Vs applied to whole model, default is sqrt(3)
%       minsta = minimum number of unique stations
%       minpha = minimum number of total phases
%           CI = confidence intervals to report (e.g. 0.95 for 95%)
%       
%
%  OUTPUTS
%
%      hyp = [lat,lon,dep] estimated hypocenter
%        T = estimated origin time (a datetime)
%    picks = [[picks], res (Tpred-Tobs), D epicentral, azimuth, take-off angle]
%    stats = a dictionary with msft,msft0...


% -- Sort input/default parameters
params   = unpack_paramsHypoAPP(params);

% -- Check that station/phase requirements are met
Np0 = size(picks,1);
Ns0 = size(unique(picks(:,1:2),'rows'),1);
if Np0 < params.minpha || Ns0 < params.minsta
    stats.Ns0   = Ns0;
    stats.Np0   = Np0;
    return
end

% -- Check that depth is within bounds
hyp( hyp(:,3) < params.minZ ,3) = params.minZ;
hyp( hyp(:,3) > params.maxZ ,3) = params.maxZ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     Get initial residuals and misfit        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tObs  = picks(:,4);
d     = distance(hyp(1),hyp(2),picks(:,1),picks(:,2),[6371, 0]);
tPred = (params.vpvs.^(picks(:,3)-1)).*RayTrace(hyp(3),d,model);
wght  = Np0*picks(:,5)/sum(picks(:,5));
res0  = (tObs-tPred)-WeightedMedian(tObs-tPred,wght);
msft0 = sum(wght.*abs(res0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%      Begin iterative grid search       %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:params.max_iter

    % -- "hypoAPPgrid" does an iteration. Store results as temporary 
    % -- variables in case misfit increases
    [hyp1,T,picks,stats] = hypoAPPgrid(hyp,T,picks,model,params);
    
    % -- Compute length of hypocenter shift
    dH = sqrt(distance(hyp(1),hyp(2),hyp1(1),hyp1(2),[6371 0])^2 ...
                   + (hyp1(3)-hyp(3))^2);     
    hyp = hyp1; 
    
    % -- If misfit did not decrease, accept step but stop iterating
    if stats.msft >= stats.msft0
        break
    end
    
    % -- Check for outlier picks
    res    = picks(:,6);
    stdres = std(res);
    jbad   = find(abs(res)>params.picktol*std(res));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%       Remove any outlier picks       %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(jbad)
        if length(jbad>1)
            [~,jbad] = max(abs(res));
        end
        picks(jbad,:) = [];
        
        % -- Check if minimum station/phase reqquirements are still met
        Np = size(picks,1);
        Ns = size(unique(picks(:,1:2),'rows'),1);
        if Np < params.minpha || Ns < params.minsta
            break 
        end
        continue
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%   If there are no outliers, check if converged    %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    else      
        
        % -- It's converged if change in hypocenter AND bounds are < tol
        if max([dH, abs(stats.rENZ-params.rENZ)]) < params.tol
            break
        else
            params.rENZ = stats.rENZ;
        end           
    end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       Add initial condidtions into stats       %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stats.msft0 = msft0;
stats.res0  = res0;
stats.Ns0   = Ns0;
stats.Np0   = Np0;
stats.iter  = iter;
        