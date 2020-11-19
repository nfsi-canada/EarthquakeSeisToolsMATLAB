% 2020-05-13
% This is the driver from a "hypoTD" test, modified from my previous 
% "synRelMT" test. There are several leftover features from the relative
% MT aspect that are not used here...I've tried to label them.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Synthetic parameters

ne = 16; % number of events
ns = 10;  % number of stations
no = 4;  % number of 'octants' with stations 
         % if e.g. no = 2, stations are restricted to upper-east quadrant
         % of focal sphere.
EvSp = 0.1;  % Event 'separation', 0 indicates truly colocated events, in km
EvDp = 20; % Mean event depth in km
EvHE = 0.5;  % ESTIMATED hypocenter error, in km
EvHS = 0;  % ESTIMATED hypocenter centroid shift, in km

jref = 1;

% Noise 0: remove phases (portions to remove [P,S])
% Noise 1: 
% Noise 2: 
% Noise 3: 
NF0 = [0.0, 0.0];
NF1 = 0;
NF2 = 0;
NF3 = 0;

model = 'ma2011';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:\Users\alex\Nextcloud\MATLAB\hypoTD\synHypoTD_pieces')
addpath('C:\Users\alex\Nextcloud\MATLAB\hypoTD\synHypoTD_functions\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Make synthetic data

% -- I put the creation of
createEventsHypoTD;
createStationsHypoTD;
createPhadirHypoTD;
sP = [];
% -- synRelMT didn't use lat/lon coordinates...need to add that here

%jE = SelectEventCombinations(hyp0,phadir,Inf,ne);
%mcc = phadir2dt(phadir,hyp0,jE);

params.vpvs = sqrt(3);
params.wCC = 0;
params.wCT = 1;
params.wsP = 0;
params.w0  = 5e-6;
params.w00 = 0;
params.wMZ = 0;
params.minZ = 0.1;
params.maxZ = 200;
params.max_iter = 1000;

[hyp,T,stats] = hypoTD(hyp0,T0,staLL,phadir,[],sP,model);

map = m_proj('UTM','lat',[min(hyp(:,1))-0.1,max(hyp(:,1))+0.1],'long',[min(hyp(:,2))-0.1,max(hyp(:,2))+0.1],'ell','wgs84','rect','on');
[exF,eyF] = m_ll2xy(hyp(:,2),hyp(:,1));
ezF = hyp(:,3);


exF = (exF-mean(exF))/1000  + mean(ex);
eyF = (eyF-mean(eyF))/1000  + mean(ey);
ezF =  ezF-mean(ezF)        + mean(ez);

plotEvStaHypoTD;

%plotStationAngles;

