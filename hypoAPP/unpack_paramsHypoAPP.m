function [R,Z,nR,nZ,minZ,maxZ,picktol,vpvs,max_iter] = unpack_paramsHypoAPP(params)
% function [R,Z,nR,nZ,minZ,maxZ,picktol,vpvs,max_iter] = unpack_paramsHypoAPP(params)
%
% 2020-06-25
% This function just unpacks optional paramaters fed to hypoAPP through the
% "params" structure, or if not provided uses the defaults. 
%
% Fields in params are:
%
%            R = maximum epicentral distance (km)
%            Z = maximum depth difference (km)
%           nR = number of horizontal grid points 
%           nZ = number of depth grid points
%                It makes sense for nR,nZ to be odd numbers so that hyp0
%                is one of the grid points
%         minZ = minimum allowable depth (km)
%         maxZ = maximum allowable depth (km)
%      picktol = cull picks with residual > picktol times the std. dev.
%         vpvs = constant Vp/Vs applied to whole model, default is sqrt(3)
%     max_iter = maximum number of iterations hypoAPP will perform

DEFparams.R  = 200;
DEFparams.Z  = 8;
DEFparams.nR = 21;
DEFparams.nZ = 9;
DEFparams.minZ = 0.1;
DEFparams.maxZ = 800;
DEFparams.picktol = 2.5;
DEFparams.vpvs = sqrt(3);
DEFparams.max_iter = 50;

if ~exist('params') 
    params = DEFparams;
else
    field = fieldnames(params);
    for ii = 1:length(field)
        DEFparams.(field{ii}) = params.(field{ii});
    end
    params = DEFparams;
end
    
R       = params.R;
Z       = params.Z;
nR      = params.nR;
nZ      = params.nZ;
minZ    = params.minZ;
maxZ    = params.maxZ;
picktol = params.picktol;
vpvs    = params.vpvs;
max_iter = params.max_iter;