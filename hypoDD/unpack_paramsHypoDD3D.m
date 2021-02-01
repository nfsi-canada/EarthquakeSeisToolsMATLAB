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
%        Niter = number of iterations
%         minZ = minimum allowable depth (km)
%         maxZ = maximum allowable depth (km)
%         minV = enforce minimum P velocity when interpolating the  
%                velocity model for ray tracing (km s-1)
%         maxV = '' maximum velocity (km s-1)
%           wP = relative weight of P data
%           wS = relative weight of S data
%          wCT = relative weight of pick-based differential times
%          wCC = relative weight of waveform-based differenctial times
%                all the weights (wP,wS,wCT,wCC) can be constant or 
%                length Niter
%         lmbd = damping factor (can be Niter x 1 vector)
%       alphCT = cutoff threshold (# of std. dev.) to cull CT equations
%       alphCC = cutoff threshold (# of std. dev.) to cull CC equations
%                   alphCT,alphCC can be Niter x 1 vectors
%          W00 = weight of zero-shift equations
%       minsta = minimum number stations with either a P or S pick
%     minstaPS = minimum  (# sta w/ P) +  (# sta w/ S)
%  NewtonSteps = max. number of Newton steps in each iteration's
%                "solveHypoDD" call
%      CGsteps = max. nubmer of conjugate-gradient steps to find
%                Newton-direction, for each Newton step
%        NobsC = min. number of unique stations to have form a continuous 
%                cluster
%        Nboot = number of bootstrap iterations (default to 0)
%          wd0 = max. inter-event distance (km) for full weight.
%                further than this comparisons are down-weighted
%                by a linear function that extend another wd1 km, to a 
%                minimum weight of 0.1
%          wd1 = linear portion of distance weighting (km)


DEFparams.vpvs     = sqrt(3);
DEFparams.Niter    = 10;
DEFparams.minZ     = 0.1;
DEFparams.maxZ     = 800;
DEFparams.minV     = 3;
DEFparams.maxV     = 10;
DEFparams.wP       = 1;
DEFparams.wS       = 1; 
DEFparams.wCT      = 1; 
DEFparams.wCC      = 1; 
DEFparams.lmbd     = 2e-4;
DEFparams.alphCT   = 5;
DEFparams.alphCC   = 5;
DEFparams.W00      = 100;
DEFparams.minsta   = 4;
DEFparams.minstaPS = 4;
DEFparams.NewtonSteps = 10;
DEFparams.CGsteps     = 20;
DEFparams.NobsC    = 6;
DEFparams.Nboot    = 0;
DEFparams.wd0      = 4;
DEFparams.wd1      = 10;

if ~exist('params') 
    params = DEFparams;
else
    field = fieldnames(params);
    for ii = 1:length(field)
        DEFparams.(field{ii}) = params.(field{ii});
    end
    params = DEFparams;
end

% -- If constant weights are given, apply them to each iteration 
if length(params.wP) == 1
    params.wP = params.wP*ones(params.Niter,1);
end
if length(params.wS) == 1
    params.wS = params.wS*ones(params.Niter,1);
end
if length(params.wCT) == 1
    params.wCT = params.wCT*ones(params.Niter,1);
end
if length(params.wCC) == 1
    params.wCC = params.wCC*ones(params.Niter,1);
end
if length(params.lmbd) == 1
    params.lmbd = params.lmbd*ones(params.Niter,1);
end
if length(params.alphCT) == 1
    params.alphCT = params.alphCT*ones(params.Niter,1);
end
if length(params.alphCC) == 1
    params.alphCC = params.alphCC*ones(params.Niter,1);
end
if length(params.wd0) == 1
    params.wd0 = params.wd0*ones(params.Niter,1);
end
if length(params.wd1) == 1
    params.wd1 = params.wd1*ones(params.Niter,1);
end

% -- Normalize relative weights of CT/CC and P/S data at each iteration
mnPS = (params.wP+params.wS)/2;
mnTC = (params.wCT+params.wCC)/2;
params.wP  = params.wP./mnPS;
params.wS  = params.wS./mnPS;
params.wCT = params.wCT./mnTC;
params.wCC = params.wCC./mnTC;

%      picktol = cull picks with residual > picktol times the std. dev.
%      steplen = initial step length to try (km)
%DEFparams.picktol = 2.5;
%DEFparams.steplen = 2;