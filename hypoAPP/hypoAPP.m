function [hyp,T,picks,stats] = hypoAPP(hyp,T,picks,model,params)
% function [hyp,T,picks,stats] = hypoAPP(hyp,T,picks,model,params)
%
% 2020-03-14
% This function estimates an earthquake hypocenter by grid-seaching around
% an initial hypocenter estimate (calling the function hypoAPPgrid). It
% repeats this process
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
%            R = initial maximum epicentral distance (km)
%            Z = initial maximum depth difference (km)
%           nR = number of horizontal grid points 
%           nZ = number of depth grid points
%                It makes sense for nR,nZ to be odd numbers so that hyp0
%                is one of the grid points
%         minZ = minimum allowable depth (km)
%         maxZ = maximum allowable depth (km)
%          tol = stop iterating when grid size less than "tol" km
%      picktol = cull picks with residual > picktol * std. dev.
%         vpvs = constant Vp/Vs applied to whole model, default is sqrt(3)
%       
%
%  OUTPUTS
%
%      hyp = [lat,lon,dep] estimated hypocenter
%        T = estimated origin time (a datetime)
%    picks = [[picks], res (Tpred-Tobs), D epicentral, azimuth, take-off angle]
%    stats = a dictionary with msft,msft0...


iter = 0;
npk0 = size(picks,1);

% -- Sort input/default parameters
[R,Z,nR,nZ,minZ,maxZ,picktol,vpvs,max_iter] = unpack_paramsHypoAPP(params);

% -- Check that depth is within bounds
hyp( hyp(:,3) < minZ ,3) = minZ;
hyp( hyp(:,3) > maxZ ,3) = maxZ;

% -- hypoAPP only takes in 1 horizontal search distance, but it adjusts
% -- E-W and N-S distances separately.
ENZ   = [R,R,Z];
CNVRG = [0 0 0];

while iter < max_iter

    hyp_prev = hyp;
  
    % -- "hypoAPPgrid" does an iteration. Store results as temporary 
    % -- variables in case misfit increases
    [hyp1,T1,picks1,stats1] = hypoAPPgrid(hyp,T,picks,model,ENZ(1),ENZ(2),ENZ(3),nR,nZ,minZ,maxZ,vpvs);
    
    % -- On the first iteration, store original misfit, and make a 
    % -- version of "stats" to ouput in case misfit already increased
    if ~iter
        msft00 = stats1.msft0;
        stats = stats1;
        stats.msft = stats.msft0;
        stats.res  = stats.res0;
    end
    
    % -- If misfit did not decrease, don't accept step and stop iterating
    if stats1.msft >= stats1.msft0
        break
    end
   
    hyp    = hyp1;
    T      = T1;
    picks  = picks1;
    stats  = stats1;
    iter   = iter+1;
    
    % -- Check for outlier picks
    res    = picks(:,6);
    stdres = std(res);
    jbad   = find(abs(res)>picktol*std(res));
    if length(jbad)
        if length(jbad>1)
            [~,jbad] = max(abs(res));
        end
        picks(jbad,:) = [];
        continue
    else    
      
        ENZt = stats.ENZ;        
        CNVRG(ENZt <  ENZ*0.99) = 0;
        CNVRG(ENZt >= ENZ*0.99) = 1;
        
        if prod(CNVRG) 
            break
        else
            ENZ = ENZt;
        end
            
    end
        
end
stats.msft0 = msft00;
stats.iter  = iter;
        