function params = unpack_paramsHypoAPP(params)
% function params = unpack_paramsHypoAPP(params)
%
% 2020-06-25
% This function just unpacks optional paramaters fed to hypoAPP through the
% "params" structure, or if not provided uses the defaults. 
%
% Fields in params are:
%
%         rENZ = initial maximum distance in [E,N,Z] directions (km)
%           NR = number of horizontal grid points 
%           NZ = number of depth grid points
%                It makes sense for nR,nZ to be odd numbers so that hyp0
%                is one of the grid points
%         minZ = minimum allowable depth (km)
%         maxZ = maximum allowable depth (km)
%      picktol = cull picks with residual > picktol times the std. dev.
%         vpvs = constant Vp/Vs applied to whole model, default is sqrt(3)
%     max_iter = maximum number of iterations hypoAPP will perform
%       minsta = minimum number of unique stations
%       minpha = minimum number of total phases
%           CI = confidence intervals to report (e.g. 0.95 for 95%)

DEFparams.rENZ     = [100,100,20];
DEFparams.NH       = 31;
DEFparams.NZ       = 27;
DEFparams.minZ     = 0.1;
DEFparams.maxZ     = 800;
DEFparams.picktol  = 2.5;
DEFparams.vpvs     = sqrt(3);
DEFparams.max_iter = 50;
DEFparams.minsta   = 4; 
DEFparams.minpha   = 6;
DEFparams.CI       = 0.90;

if ~exist('params') 
    params = DEFparams;
else
    field = fieldnames(params);
    for ii = 1:length(field)
        DEFparams.(field{ii}) = params.(field{ii});
    end
    params = DEFparams;
end
    
% -- Enure NR, NZ are odd numbers
params.NH = params.NH + 1-mod(params.NH,2);
params.NZ = params.NZ + 1-mod(params.NZ,2);


