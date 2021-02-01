function err = bootstrapHypoTD(G,dCT,dCC,dST,dSC,dZ,z,params);
% function err = bootstrapHypoTD(G,dCT,dCC,dST,dSC,dZ,z,params);
%
% 2020-01-14
% This function takes in a full hypoTD linear system, which should be 
% computed from the final hypocenter estimates, and bootstraps to estimate
% error in the hypocenters. Each bootstrap iteration resamples the
% residuals (dCT,dCC) with replacement.
%
% This was modified from "bootstrapHypoDD.m"
%
%
% INPUTS
%
%     G == hypoTD matrix formed with:
%         GCT == pick-based TD matrix (formed with uniform weight)
%         GCC == waveform-based TD matrix (formed with uniform weight)
%         GST == pick-based station DD matrix (formed with uniform weight)
%         GSC == waveform station DD matrix (formed with uniform weight)
%          GZ == zero-shift constraint matrix
%     dCT == pick-based DD residuals (formed with uniform weight)
%     dCC == waveform-based DD residuals (formed with uniform weight)
%     dST == pick-based station DD residuals (uniform weight)
%     dSC == waveform-based station DD residuals (uniform weight)
%            dCT,dCC,dST,dCC are randomly resampled (independently)
%      dZ == zero-shift constraint RHS vector (zeros)
%       z == earthquake depths of final
%  params ==  provided to hypoDD) the parameters we need here are
%                Nboot == number of boot
%                NewtonSteps == max. number of steps for each solution
%                CGsteps == number of conjugate-gradient steps for each
%                           Newton step
%                prct = percentile (e.g. 90) to report
%
% OUTPUTS 
%
%    err == Ne x 3 matrix [errX,errY,errZ] of 95% confidence errors
%           spatial errors in km, time error in s
%%%%%%%%%
%  errV == 3 x 3 X Ne tensor of principal error axes (true to scale)



Nct = length(dCT);
Ncc = length(dCC);
Nst = length(dST);
Nsc = length(dSC);
Ne  = size(G,2)/3;
m   = zeros(params.Nboot,3*Ne);

for ii = 1:params.Nboot

   dCTi = dCT(randi(max(Nct,1),[Nct,1]));
   dCCi = dCC(randi(max(Ncc,1),[Ncc,1]));
   dSTi = dST(randi(max(Nst,1),[Nst,1]));
   dSCi = dSC(randi(max(Nsc,1),[Nsc,1]));
   d    = [dCTi; dCCi; dSTi; dSCi; dZ];

   % -- Don't worry about depth constraints here, they might make
   % -- errors look artificially small
   m(ii,:) = solveHypoTD(G,d,z,params.NewtonSteps,params.CGsteps,-Inf,Inf);
end

% -- Take 90th percentile of [errX,errY,errZ,errT] (km or s)
err = prctile(abs(m),90)';
err = reshape(err,Ne,3);


