function err = bootstrapHypoDD(G,dCT,dCC,dZ,z,params);
% function err = bootstrapHypoDD(G,dCT,dCC,dZ,z,params);
%
% 2020-01-14
% This function takes in a full hypoDD linear system, which should be 
% computed from the final hypocenter estimates, and bootstraps to estimate
% error in the hypocenters. Each bootstrap iteration resamples the
% residuals (dCT,dCC) with replacement.
%
%
% INPUTS
%
%     G == hypoDD matrix formed with:
%         GCT == pick-based DD matrix (formed with uniform weight)
%         GCC == waveform-based DD matrix (formed with uniform weight)
%          GZ == zero-shift constraint matrix
%     dCT == pick-based DD residuals (formed with uniform weight)
%     dCC == waveform-based DD residuals (formed with uniform weight)
%            dCT,dCC are what are randomly resampled (independently)
%      dZ == zero-shift constraint RHS vector (zeros)
%       z == earthquake depths of final
%  params ==  provided to hypoDD) the parameters we need here are
%                Nboot == number of boot
%                NewtonSteps == max. number of steps for each solution
%                CGsteps == number of conjugate-gradient steps for each
%                           Newton step
%
% OUTPUTS 
%
%    err == Ne x 4 matrix [errX,errY,errZ,errT] of 95% confidence errors
%           spatial errors in km, time error in s
%%%%%%%%%
%  errV == 3 x 3 X Ne tensor of principal error axes (true to scale)



Nct = length(dCT);
Ncc = length(dCC);
Ne  = size(G,2)/4;
m   = zeros(params.Nboot,4*Ne);

for ii = 1:params.Nboot

   dCTi = dCT(randi(max(Nct,1),[Nct,1]));
   dCCi = dCC(randi(max(Ncc,1),[Ncc,1]));
   d    = [dCTi; dCCi; dZ];

   % -- Don't worry about depth constraints here, they might make
   % -- errors look artificially small
   m(ii,:) = solveHypoDD(G,d,z,params.NewtonSteps,params.CGsteps,-Inf,Inf);
end

% -- Take 95th percentile of [errX,errY,errZ,errT] (km or s)
err = prctile(abs(m),95)';
err = reshape(err,Ne,4);

% -- It would be good to record error ellipse info...
% -- How to convert SVD scale to something physically intuitive???
% errV = zeros(3,3,Ne);
% for ie = 1:Ne
%     mE = m(:,ie+Ne*(0:2));
%     [~,S,V] = svd(mE,0);
%     errV(:,:,ie) = S*V;
% end
% 


   %d    = [dCTi; dCCi; dZ; dD];

% -- Double check weight of damping matrix...it needs to be minimal
% -- for error analysis to be at all valid. Make it just high enough for 
% -- Numerical stability.
%j1D  = size(G,1)+(-length(dD)+1:0)';
%GD   =  G(j1D,:);
%G(j1D,:) = GD*min(100*eps/min(GD(find(GD))));