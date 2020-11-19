% 2020-05-06
% This piece of hypoTD builds a matrix to impose sP depth constraints


if isempty(sP)
    AsP = zeros(0,3*Ne);
    bsP = zeros(0,1);
else
   % dES
   
   % -- Build the simple matrix of ones
   AsP = sparse(1:NsP,2*Ne+sP(:,1),ones(NsP,1),NsP,3*Ne);
   bsP = zeros(NsP,1);
   
   % -- This method might be able to be more efficient...I'm making grids
   %    that are very similar for one event to the next...
   for ii = 1:NsP
       
       % -- Epicentral distance (km) and sP time (seconds after P)
       Xsp = dES(sP(ii,1),sP(ii,2));
       Tsp = sP(ii,3);
           
       % -- Need times to a range of "bounce points", with the maximum
       % -- at the current estiamte for epicentral distance (EQ to station)
       % -- D,Z are grids of bounce-point, EQ depth (x,z)
       xrngSP = linspace(1e-6,Xsp,200)';
       
       [Zt,Dt]  = meshgrid(zrng,xrngSP);
       
       % -- Compute upward P-time, upward S-time, and direct P-time
       % -- (as a function of bounce-point and depth)
       tp1 = interp2(Zrng,Drng,tp,Zt,Dt);
       ts1 = tp1*vpvs;
       tp0 = tp1(end,:);
       
       % -- Post-Bounce Pwave -time as a function of bounce-point
       tp2 = interp1(drng,tBounce,flipud(xrngSP));
       tp2 = repmat(tp2,1,length(xrngSP));
       
       % -- Total sP (realtive to direct P) time as a function of 
       % -- bounce-point and depth. Minimize w.r.t bounce point.
       gtsP = ts1 + tp2 - tp0;
       gstP(isnan(gtsP)) = 0;
       tsPz = min(gtsP,[],1);
       
       % -- Estimate depth from sP time
       try
       bsP(ii) = interp1(tsPz,zrng,Tsp);
       catch
           keyboard
       end
        
   end    
 
end