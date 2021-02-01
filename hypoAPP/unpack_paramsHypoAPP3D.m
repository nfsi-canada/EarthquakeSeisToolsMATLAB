function params = unpack_paramsHypoAPP3D(params)
% function params = unpack_paramsHypoAPP3D(params)
%
% 2020-06-25
% This function just unpacks optional paramaters fed to hypoAPP3D through 
% the "params" structure, or if not provided uses the defaults. 
%
% Fields in params are:
%
%         vpvs = constant Vp/Vs applied to whole model, default is sqrt(3)
%          tol = minimum step length to worry about (km)
%      picktol = cull picks with res. > picktol times med. abs. deviation
%     max_iter = maximum number of iterations hypoAPP will perform
%         minZ = minimum allowable depth (km)
%         maxZ = maximum allowable depth (km)
%       minpha = minimum number of picks (P+S) 
%        Lcube = side length (km) of cube to perform 'grid-search' over
%     Nsteptol = the programs stops iterating if the distance between
%                the current hypocenter and the one from Nsteptol 
%                iterations ago are with "tol" (km)
%           CI = decimal confidence intervals to return, e.g. 0.95 
%                to return 95% confidence intervals
%         minV = enforce minimum P velocity when interpolating the  
%                velocity model for ray tracing (km s-1)
%         maxV = '' maximum velocity (km s-1)
%         sclF = If the cube being divided has side length <= sclF (km)
%                hypoAPP3D will only call RayTrace3D once, and just 
%                perturb the ray-path 


DEFparams.vpvs     = sqrt(3);
DEFparams.tol      = 0.1;
DEFparams.picktol  = 6;
DEFparams.max_iter = 50;
DEFparams.minZ     = 0.1;
DEFparams.maxZ     = 800;
DEFparams.minpha   = 6;
DEFparams.Lcube    = 50;
DEFparams.NstepTol = 10;
DEFparams.CI       = 0.90;
DEFparams.minV     = 3;
DEFparams.maxV     = 10;
DEFparams.sclF     = 6; 

if ~exist('params') 
    params = DEFparams;
else
    field = fieldnames(params);
    for ii = 1:length(field)
        DEFparams.(field{ii}) = params.(field{ii});
    end
    params = DEFparams;
end
    
%      picktol = cull picks with residual > picktol times the std. dev.
%      steplen = initial step length to try (km)
%DEFparams.picktol = 2.5;
%DEFparams.steplen = 2;