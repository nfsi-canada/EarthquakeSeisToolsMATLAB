% 2021-01-15
% This is the driver for a triple-difference hypocenter inversion test,
% modified from "synHypoDD" on 2021-01-15

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Synthetic parameters

Ne = 24; % number of events
Ns = 6;  % number of stations
No = 4;  % number of 'octants' with stations 
         % if e.g. no = 2, stations are restricted to upper-east quadrant
         % of focal sphere.
NsP = 4; % number of sP depth constraints
         
scl  = 3;    % Event 'separation', 0 indicates truly colocated events,
             % standard deviation of EQ coordintates in km
errR  = 2.0; % gausian hypocenter error S.D. in km
errS  = 0;   % mean hypocenter error centroid shift, in km
errSP = 0.0; % gaussian of time error added to sP-P times (s)

map  = m_proj('lambert', 'lat',[42,68],'long',[-180,-140],'ell','wgs84','rect','on');
latC = 55;
lonC = -157;
depC = 40;
vpvs = 1.768;

% -- "Ma2011" velocity model
model = [  0.0  6.25
           8.0  6.50
          17.0  6.60
          24.0  6.70
          30.0  7.10
          36.0  7.20
          37.0  8.00
          60.0  8.05
         185.0  8.10
         220.0  8.55
         600.0 10.15];
     
% Noise 0: remove phases (portions to remove [P,S])
% Noise 1: standard devation of noise to picks (seconds)
NF0 = [0.0, 0.0];
NF1 = 0.00;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Make synthetic data

% -- 
[hyp0,hyp1] = createEventsHypoDD(Ne,latC,lonC,depC,scl,errR,errS);
T0 = datetime(2021,01,01)+days(0:Ne-1)';

stadir = createStationsHypoDD(Ns,No,model,latC,lonC,depC);
[phadir0,phadir1,T1] = createPhadirHypoDD(hyp0,hyp1,T0,stadir,model,vpvs,NF0,NF1);

% -- Create sP data
params.msp = create_sP_hypoDD(NsP,errSP,hyp0,hyp1,stadir,phadir0,model,vpvs);

% -- Form double-difference times
% -- ph2dt(hyp,phadir,maxD,maxN)
params.mct = ph2dt(hyp1,phadir1,20,12);
params.mcc = [];

% -- Form triple-differences
params.mcttd = TripleDifferenceFull(params.mct);
params.mcctd = TripleDifferenceFull(params.mcc);

% -- Form station-pair double-difference times
params.mst = ph2dtSTA(phadir1);
params.msc = [];


% -- Set parameters
params.vpvs     = vpvs;
params.Niter    = 10;
params.minZ     = 0.1;
params.maxZ     = 800;
params.minV     = 3;
params.maxV     = 10;
params.wP       = 1;
params.wS       = 1; 
params.wCT      = 1; 
params.wCC      = 1; 
params.wST      = 1; 
params.wSC      = 1;
params.lmbd     = 2e-4;
params.alphCT   = 6;
params.alphCC   = 6;
params.W00      = 0.1;
params.minsta   = 4;
params.minstaPS = 4;
params.Nboot    = 200;
params.NobsC    = 3;

% -- Actually run hypoDD!!!
[hyp,T,stats] = hypoTDX(hyp1,T1,stadir,phadir1,model,params);

jr = find(~isnan(hyp(:,1)));
[ex,ey]   = m_ll2xy(hyp(jr,2),hyp(jr,1));
[ex0,ey0] = m_ll2xy(hyp0(jr,2),hyp0(jr,1));
[ex1,ey1] = m_ll2xy(hyp1(jr,2),hyp1(jr,1));
ez  = hyp(jr,3);
ez0 = hyp0(jr,3);
ez1 = hyp1(jr,3);

ex0 = ex0/1000;
ey0 = ey0/1000;
ex1 = ex1/1000;
ey1 = ey1/1000;
ex  = (ex-mean(ex))/1000  + mean(ex0);
ey  = (ey-mean(ey))/1000  + mean(ey0);
ez  =  ez-mean(ez)        + mean(ez0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Plot events

figure(30289)
clf
hold on

plot3(ex0,ey0,ez0,'ko','MarkerFaceColor',[0.5,0.5,0.5])
plot3(ex,ey,ez,'ko','MarkerFaceColor','b')
plot3([ex0';ex'],[ey0';ey'],[ez0';ez'],'r--')

set(gca,'ZDir','reverse')
grid on
axis equal

xlim([min([ex0;ex])-0.5, max([ex0;ex])+0.5])
ylim([min([ey0;ey])-0.5, max([ey0;ey])+0.5])
zlim([min([ez0;ez])-0.5, max([ez0;ez])+0.5])

xlabel('E (km)')
ylabel('N (km)')
zlabel('Depth (km)')

errX = sqrt((ex0-ex).^2 + (ey0-ey).^2 + (ez0-ez).^2 ); 
errX1 = sqrt((ex0-ex1).^2 + (ey0-ey1).^2 + (ez0-ez1).^2 ); 
title(['errX,errX1 = [',num2str(mean(errX),'%5.3f'),', ',num2str(mean(errX1),'%5.3f'),'] km'])




