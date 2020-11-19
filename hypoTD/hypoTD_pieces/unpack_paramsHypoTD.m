% 2020-05-06
% This piece of hypoTD just handles some parameters that may or may not 
% have been specified in the structure "params"
%
%

DEFparams.wCC = 2e1;
DEFparams.wCT = 5e0;
DEFparams.wsP = 5e2;
DEFparams.w0  = 1e-1;
DEFparams.max_iter = 15;
DEFparams.vpvs = sqrt(3);
DEFparams.NsgCC = 3;
DEFparams.NsgCT = 3;
DEFparams.NsgSP = 3;
DEFparams.minZ  = 1e-3;
DEFparams.maxZ  = 350;
DEFparams.mxDCC = 10; 
DEFparams.mxDCT = 10;
DEFparams.mxNB  = 10;  % max. number of neighbors 
DEFparams.mxTD  = 12;

if ~exist('params') 
    params = DEFparams;
else
    field = fieldnames(params);
    for ii = 1:length(field)
        DEFparams.(field{ii}) = params.(field{ii});
    end
    params = DEFparams;
end
    
wCC = params.wCC;
wCT = params.wCT;
wsP = params.wsP;
w0  = params.w0;
max_iter = params.max_iter;
vpvs = params.vpvs;
NsgCC = params.NsgCC;
NsgCT = params.NsgCT;
NsgSP = params.NsgSP;
minZ  = params.minZ;
maxZ  = params.maxZ;
mxDCC = params.mxDCC;
mxDCT = params.mxDCT;
mxNB  = params.mxNB;
mxTD  = params.mxTD;

if length(mxDCC)==1
    mxDCC = repmat(mxDCC,max_iter,1);
end
if length(mxDCT)==1
    mxDCT = repmat(mxDCT,max_iter,1);
end
if isempty(mcc)
    mcc = zeros(0,6);
end
if isempty(sP)
    sP = zeros(0,3);
end