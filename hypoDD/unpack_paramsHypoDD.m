function params = unpack_paramsHypoDD(params)
% function params = unpack_paramsHypoDD(params)
%
% 2020-06-25
% This function just unpacks optional paramaters fed to hypoAPP3D through 
% the "params" structure, or if not provided uses the defaults. 
%
% Fields in params are:
%
%         vpvs = constant Vp/Vs applied to whole model, default is sqrt(3)
%        mcttd = catalog differential times ((t1a-t1b)-(t2a-t2b))_obs
%                [EVa,EVb,STA1,STA2,PHA,DDT,(CC)]
%        mcctd = waveform differential times ((t1a-t1b)-(t2a-t2b))_obs
%                [EVa,EVb,STA1,STA2,PHA,DDT,CC]
%          msp = [EV,STA,T] sP-P times
%          mst = catalog station-pair differential times (t1-t2)_obs
%                [EV,STA1,STA2,PHA,DDT,(CC)]
%          msc = waveform station-pair differential times (t1-t2)_obs
%                [EV,STA1,STA2,PHA,DDT,(CC)]
%        mcctd = waveform differential times ((t1a-t1b)-(t2a-t2b))_obs
%                [EVa,EVb,STA1,STA2,PHA,DDT,CC]
%        Niter = number of iterations
%         minZ = minimum allowable depth (km)
%         maxZ = maximum allowable depth (km)
%           wP = relative weight of P data
%           wS = relative weight of S data
%          wCT = relative weight of pick-based differential times
%          wCC = relative weight of waveform-based differenctial times
%          wST = relative weight of pick-based station differential times
%          wSC = relative weight of waveform station differential times
%                all the weights (wP,wS,wCT,wCC,wST,wSC) can be constant or 
%                length Niter
%         lmbd = damping factor (can be Niter x 1 vector)
%       alphCT = cutoff threshold (# of std. dev.) to cull CT equations
%       alphCC = cutoff threshold (# of std. dev.) to cull CC equations
%       alphST = cutoff threshold (# of std. dev.) to cull ST equations
%       alphSC = cutoff threshold (# of std. dev.) to cull SC equations
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
%         dmax = max. inter-event distance (km) 
%                cubic weighting based on distance is used.
%          wsP = only relavant for the special "sP" versions, it weights
%                the sP depth constraints
%        prctB = percentile error to report from bootstrapping 
%                default is 95;

DEFparams.vpvs     = sqrt(3);
DEFparams.mcttd    = zeros(0,7);
DEFparams.mcctd    = zeros(0,7);
DEFparams.msp      = zeros(0,3);
DEFparams.mst      = zeros(0,6);
DEFparams.msc      = zeros(0,6);
DEFparams.Niter    = 10;
DEFparams.minZ     = 0.1;
DEFparams.maxZ     = 800;
DEFparams.wP       = 1;
DEFparams.wS       = 1; 
DEFparams.wCT      = 1; 
DEFparams.wCC      = 1; 
DEFparams.wST      = 1; 
DEFparams.wSC      = 1; 
DEFparams.lmbd     = 2e-4;
DEFparams.alphCT   = 5;
DEFparams.alphCC   = 5;
DEFparams.alphST   = 5;
DEFparams.alphSC   = 5;
DEFparams.W00      = 100;
DEFparams.minsta   = 4;
DEFparams.minstaPS = 4;
DEFparams.NewtonSteps = 10;
DEFparams.CGsteps     = 20;
DEFparams.NobsC    = 6;
DEFparams.Nboot    = 200;
DEFparams.GXmax    = 10000;
DEFparams.GZmax    = 500;
DEFparams.dmax     = 12;
DEFparams.wSP      = 100;
DEFparams.prctB    = 95;

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
if length(params.wST) == 1
    params.wST = params.wST*ones(params.Niter,1);
end
if length(params.wSC) == 1
    params.wSC = params.wSC*ones(params.Niter,1);
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
if length(params.alphST) == 1
    params.alphST = params.alphST*ones(params.Niter,1);
end
if length(params.alphSC) == 1
    params.alphSC = params.alphSC*ones(params.Niter,1);
end
if length(params.dmax) == 1
    params.dmax = params.dmax*ones(params.Niter,1);
end
if length(params.wSP) == 1
    params.wSP = params.wSP*ones(params.Niter,1);
end


% -- Make sure main data arrays weren't given as []
if ~length(params.mst)
    params.mct = zeros(0,6);
end
if ~length(params.msc)
    params.mcc = zeros(0,6);
end
if ~length(params.mst)
    params.mst = zeros(0,6);
end
if ~length(params.msc)
    params.msc = zeros(0,6);
end
if ~length(params.mcttd)
    params.mcttd = zeros(0,7);
end
if ~length(params.mcctd)
    params.mcctd = zeros(0,7);
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