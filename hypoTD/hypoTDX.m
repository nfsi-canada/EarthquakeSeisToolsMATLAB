function [hyp1,T1,stats] = hypoTDX(hyp0,T0,sta,phadir,model,params)
% function [hyp1,T1,stats] = hypoTDX(hyp0,T0,sta,phadir,model,params)
%
% 2020-01-23
% This code was modifed from "hypoDD.m", which mimics the program hypoDD
% by Felix Waldhauser. It uses triple-differences rather than
% than double-differences, removing origin times from the system (they are
% computed at the end).
%
% This is the special "X" version of "hypoTD" that can take in sP and 
% station double-difference constraints, but it just takes in these data 
% in "params".
%
% This should make hypoTD and hypoTDsP obsolete...
%
%   INPUTS 
%
%     hyp0 = [lat0,lon0,dep0] initial hypocenters (NEx3)
%       T0 = intial origin times (NEx1 datetime)
%      sta = [latS,lonS,depS(km)]
%   phadir = [EV,STA,PHA,T]
%    model = N_layer x 2 matrix [Depth,Velocity]
%   params = structure of optional paramters
%
%
%  OUTPUTS
%
%     hyp1 = [lat,lon,dep] estimated hypocenter (NEx3)
%       T1 = estimated origin times (NEx1 datetime)
%    stats = a dictionary with various output statistics


% -- Sort input/default parameters 
params = unpack_paramsHypoDD(params);

mcttd  = params.mcttd;
mcctd  = params.mcctd;
mst    = params.mst;
msc    = params.msc;
msp    = params.msp;
params = rmfield(params,{'mcttd','mcctd','mst','msc','msp'});

Ne  = size(hyp0,1);
Ns  = size(sta,1);
Np  = size(phadir,1);
Nct = size(mcttd,1);
Ncc = size(mcctd,1);
Nsp = size(msp,1);
Nst = size(mst,1);
Nsc = size(msc,1);


% -- MAD for gaussian noise (sigma_MAD in W2000 eq. 15)
MAD = 0.67449;

% -- Initialize model vector with origin times just set to zero
% -- Initialize vector of relocated events as full list
jE = (1:Ne)';

% -- Prep interpolants of travel-time, take-off angle, source-velocity
[Ft,Fg,Fv] = prep1DgridHypoDD(model,params.GXmax,params.GZmax);

% -- Prep interpolant of depth given distance and sP time
FZsp = prepInterpolant_sP(model,params.vpvs,params.GXmax,params.GZmax);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Set up map. NOTE: REQUIRES M_MAP PACKAGE    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- M_Map is a mapping packege for Matlab/Octave:
% -- Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB",
% -- https://www.eoas.ubc.ca/~rich/map.html 
% -- This code may not be compatible with lastest version of M_Map

lat0 = min([hyp0(:,1); sta(:,1)])-1;
lat1 = max([hyp0(:,1); sta(:,1)])+1;
lon0 = min([hyp0(:,2); sta(:,2)])-1;
lon1 = max([hyp0(:,2); sta(:,2)])+1;
map  =  m_proj('lambert', 'lat',[lat0,lat1],'long',[lon0,lon1],'ell','wgs84','rect','on');

% -- EQ and Station coordinates in map units (km)
[ex,ey] = m_ll2xy(hyp0(:,2),hyp0(:,1));
[sx,sy] = m_ll2xy(sta(:,2),sta(:,1));
H = [ex/1000,ey/1000,hyp0(:,3)];
S = [sx/1000,sy/1000,sta(:,3)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    Intialize structure to record stats    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stats.Ne0 = Ne;
stats.Ns0 = Ns;
stats.Np0 = Np;
stats.Nct0 = Nct;
stats.Ncc0 = Ncc;
stats.Nst0 = Nst;
stats.Nsc0 = Nsc;
stats.Nsp0 = Nsp;


stats.Nclst = zeros(params.Niter,1);
stats.Ne    = zeros(params.Niter,1);
stats.Nct   = zeros(params.Niter,1);
stats.Ncc   = zeros(params.Niter,1);
stats.Nst   = zeros(params.Niter,1);
stats.Nsc   = zeros(params.Niter,1);
stats.Nsp   = zeros(params.Niter,1);
stats.DH    = zeros(params.Niter,1);
stats.DZ    = zeros(params.Niter,1);
stats.DT    = zeros(params.Niter,1);
stats.OS    = zeros(params.Niter,3);
stats.prctM = zeros(params.Niter,1);
stats.iter  = zeros(params.Niter,1);
stats.jE    = cell(params.Niter,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    Print some initial stats    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['%%%%%%%%%%%%%%%   hypoTD starting ',datestr(datetime('now')),'   %%%%%%%%%%%%%%%'])
disp(['# earthquakes: ', num2str(Ne)])
disp(['# stations: ', num2str(Ns)])
disp(['# pick-based TDs: ',num2str(Nct)])
disp(['# waveform-based TDs: ',num2str(Ncc)])
disp(['# pick-based station DDs: ',num2str(Nst)])
disp(['# waveform-based station DDs: ',num2str(Nsc)])
disp(['# sP depth contraints: ',num2str(Nsp)])
disp(' ')
disp('Iter     Ne   Nclst      Nct     Ncc      Nst     Nsc    DH(km)  DZ(km)  %d_msft')
disp('--------------------------------------------------------------------')

%fprintf('% 3d  % 5d   % 7d  % 7d  % 6.3f  % 6.3f  % 6.3f  % 5.2f\n',  ... 
%             iter,Ne,Nct,Ncc,sts.mnDH,sts.mnDV,sts.mnDT,sts.prctM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%      Begin DD iterations      %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:params.Niter
  
    % -- Cut-off thresholds of residuals for this iteration, plus
    % -- distance (km) for weight = 0.1, damping parameter
    aCT  = params.alphCT(iter);
    aCC  = params.alphCC(iter);
    aST  = params.alphST(iter);
    aSC  = params.alphSC(iter);
    dmax = params.dmax(iter); 
    lmbd = params.lmbd(iter);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%      Cull and Cluster Events      %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % -- Remove equations for events too far apart...
    % -- NOT ACTUALLY REMOVING EQUATIONS RIGHT NOW!!
    % -- JUST WEIGHTING DISTANT ONES TO 0.1 ...
    %mcttd = checkDistanceHypoDD(mcttd,H,2*(wd0+wd1));
    %mcctd = checkDistanceHypoDD(mcctd,H,2*(wd0+wd1));
    
    % -- Cluster events, ensuring continuous connection
    % -- Remove comparison for events in different clusters?
    % -- Unlike hypoDD, I do this separately for each iteration
    [mcttd,mcctd,jE]   = checkMinStaHypoTD(mcttd,mcctd,jE,params.minsta,params.minstaPS);
    [C,mcttd,mcctd,jE] = clusterHypoTD(mcttd,mcctd,jE,params.NobsC);
    
    Ne  = length(jE);
    Nct = size(mcttd,1);
    Ncc = size(mcctd,1);
    
    % -- Check if we've eliminated all events!!
    if Ne < 2
        break
    end
    
    mst(~ismember(mst(:,1),jE),:) = [];
    msc(~ismember(msc(:,1),jE),:) = [];
    msp = msp(ismember(msp(:,1),jE),:);
    Nst = size(mst,1);
    Nsc = size(msc,1);
    Nsp = size(msp,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%       Do Ray Tracing      %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % -- Assign "phase numbers" to each unique event/station combo
    uES = unique([mcttd(:,[1 3]); mcttd(:,[2 3]); 
                  mcttd(:,[1 4]); mcttd(:,[2 4]); 
                  mcctd(:,[1 3]); mcctd(:,[2 3]); 
                  mcctd(:,[1 4]); mcctd(:,[2 4]);
                  mst(:,[1 2]); mst(:,[1 3]); 
                  msc(:,[1 2]); msc(:,[1 3])],'rows');   
                 
    % -- Do Ray-Tracing
    [t,g,vto] = RayTrace1DhypoDD(H(uES(:,1),:),S(uES(:,2),:),Ft,Fg,Fv);
 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Sort data prior to taking step %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,evsCT] = ismember(mcttd(:,1:2),jE);
    [~,evsCC] = ismember(mcctd(:,1:2),jE);
    [~,evsST] = ismember(mst(:,1),jE);
    [~,evsSC] = ismember(msc(:,1),jE);
    [~,evsSP] = ismember(msp(:,1),jE);

    % -- Index is phase number
    [~,jptA1] =  ismember(mcttd(:,[1,3]),uES,'rows');  
    [~,jptB1] =  ismember(mcttd(:,[2,3]),uES,'rows');
    [~,jptA2] =  ismember(mcttd(:,[1,4]),uES,'rows');  
    [~,jptB2] =  ismember(mcttd(:,[2,4]),uES,'rows');
    [~,jpcA1] =  ismember(mcctd(:,[1,3]),uES,'rows');
    [~,jpcB1] =  ismember(mcctd(:,[2,3]),uES,'rows');
    [~,jpcA2] =  ismember(mcctd(:,[1,4]),uES,'rows');
    [~,jpcB2] =  ismember(mcctd(:,[2,4]),uES,'rows');
    [~,jpst1] =  ismember(  mst(:,[1,2]),uES,'rows');
    [~,jpst2] =  ismember(  mst(:,[1,3]),uES,'rows');
    [~,jpsc1] =  ismember(  msc(:,[1,2]),uES,'rows');
    [~,jpsc2] =  ismember(  msc(:,[1,3]),uES,'rows');

    % -- Observed differential times
    dtCT = mcttd(:,6);
    dtCC = mcctd(:,6);
    dtST = mst(:,5);
    dtSC = msc(:,5);

    % -- Initialize weights from input (reweight P vs. S, CC vs. CT later)
    wCT  = mcttd(:,7);
    wCC  = mcctd(:,7);
    wST  = mst(:,6);
    wSC  = msc(:,6);
    
    % -- Indices of P, S equations, 
    jPct = find(mcttd(:,5)==1);
    jSct = find(mcttd(:,5)==2);
    jPcc = find(mcctd(:,5)==1);
    jScc = find(mcctd(:,5)==2);
    jPst = find(mst(:,4)==1);
    jSst = find(mst(:,4)==2);
    jPsc = find(msc(:,4)==1);
    jSsc = find(msc(:,4)==2);

    % -- Vector of Vp/Vs correction (for S)
    vfCT = ones(Nct,1);
    vfCC = ones(Ncc,1);
    vfST = ones(Nst,1);
    vfSC = ones(Nsc,1);
    vfCT(jSct) = params.vpvs;
    vfCC(jScc) = params.vpvs;
    vfST(jSst) = params.vpvs;
    vfSC(jSsc) = params.vpvs;

    % -- Distance between events, and corresponding weights
    wdCT = distanceWeightingHypoDD(mcttd,H,dmax);
    wdCC = distanceWeightingHypoDD(mcctd,H,dmax);
    
    % -- Relative weighting of CT vs. CC  and P vs. S data
    WCT = wCT.*wdCT*params.wCT(iter);
    WCC = wCC.*wdCC*params.wCC(iter);
    WST = wST.*params.wST(iter);
    WSC = wSC.*params.wSC(iter);
    WCT(jPct) = WCT(jPct)*params.wP(iter);
    WCT(jSct) = WCT(jSct)*params.wS(iter);
    WCC(jPcc) = WCC(jPcc)*params.wP(iter);
    WCC(jScc) = WCC(jScc)*params.wS(iter);
    WST(jPst) = WST(jPst)*params.wP(iter);
    WST(jSst) = WST(jSst)*params.wS(iter);
    WSC(jPsc) = WSC(jPsc)*params.wP(iter);
    WSC(jSsc) = WSC(jSsc)*params.wS(iter);
    
    % -- Index unique ray-paths to get take-off angles
    gCT = [g(jptA1,:) g(jptB1,:) g(jptA2,:) g(jptB2,:)];
    gCC = [g(jpcA1,:) g(jpcB1,:) g(jpcA2,:) g(jpcB2,:)];
    gST = [g(jpst1,:) g(jpst2,:)];
    gSC = [g(jpsc1,:) g(jpsc2,:)];
    
    % -- Index unique ray-paths again to get velocities at source
    % -- Divide by Vp/Vs for S-data
    vCT = [vto(jptA1)./vfCT vto(jptB1)./vfCT];
    vCC = [vto(jpcA1)./vfCC vto(jpcB1)./vfCC];
    vST = vto(jpst1)./vfST;
    vSC = vto(jpsc1)./vfSC;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  Actually build linear system and solve %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % -- Form core DD matrices and RHS vectors 
    GCT = buildHypoTDmatrix(evsCT,gCT,vCT,WCT,Ne);
    GCC = buildHypoTDmatrix(evsCC,gCC,vCC,WCC,Ne);
    dCT = WCT.*(dtCT-vfCT.*( (t(jptA1)-t(jptB1))-(t(jptA2)-t(jptB2)) ));
    dCC = WCC.*(dtCC-vfCC.*( (t(jpcA1)-t(jpcB1))-(t(jpcA2)-t(jpcB2)) ));
    
    % -- Form station DD matrices and RHS vectors
    GST = buildStationDDmatrix(evsST,gST,vST,WST,Ne,1);
    GSC = buildStationDDmatrix(evsSC,gSC,vSC,WSC,Ne,1);
    dST = WST.*(dtST-vfST.*(t(jpst1)-t(jpst2)));
    dSC = WSC.*(dtSC-vfSC.*(t(jpsc1)-t(jpsc2)));
    
    % -- Form sP depth constraint matrix
    [GsP,dsP] = sPmatrixHypoTD(msp,H,S,FZsp,params.wSP(iter),jE);
    
    % -- Form zero-shift and damping matrices and RHS vectors
    [GZ,dZ] = ZeroShiftHypoTD(C,params.W00,evsSP);
    [GD,dD] = DampingMatrixHypoTD(C,stats.Ne0,lmbd);
    
    
    % -- Combine portions of systems and solve with Newton's method
    G = [GCT; GCC; GST; GSC; GsP; GZ; GD];   
    d = [dCT; dCC; dST; dSC; dsP; dZ; dD];
    
    [m,s] = solveHypoTD(G,d,H(jE,3),params.NewtonSteps,params.CGsteps, ... 
                                    params.minZ,params.maxZ);
   
    fprintf(' % 3d   % 5d   % 3d   % 7d  % 6d  % 7d  % 6d   % 6.3f  % 6.3f   % 5.2f\n',  ... 
             iter,Ne,max(C),Nct,Ncc,Nst,Nsc,s.DH,s.DZ,s.prctM)                           

    % -- Update hypocenter array with computed optimal shift
    H(jE,1) = H(jE,1) + m(1:Ne);
    H(jE,2) = H(jE,2) + m(Ne+[1:Ne]);
    H(jE,3) = H(jE,3) + m(2*Ne+[1:Ne]);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Check for outlier equations %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iter < params.Niter
        % -- Indices, Weights for CT,CC,...
        jct = (1:Nct)';
        jcc = Nct+(1:Ncc)';
        jst = Nct+Ncc+(1:Nst)';
        jsc = Nct+Ncc++Nst+(1:Nsc)';
        
        % -- Reweight and cull equations based on residual (eq. 15 of W2000)
        drCT  = G(jct,:)*m-d(jct);
        drCC  = G(jcc,:)*m-d(jcc);
        drST  = G(jst,:)*m-d(jst);
        drSC  = G(jsc,:)*m-d(jsc);
        MADCT = median(abs(drCT-median(drCT)));
        MADCC = median(abs(drCC-median(drCC)));
        MADST = median(abs(drST-median(drST)));
        MADSC = median(abs(drSC-median(drSC)));
        
        mcttd(:,7) = max(0,1-((MAD*drCT)/(aCT*MADCT)).^2).^2;
        mcctd(:,7) = max(0,1-((MAD*drCC)/(aCC*MADCC)).^2).^2;
        mst(:,6)   = max(0,1-((MAD*drST)/(aST*MADST)).^2).^2;
        msc(:,6)   = max(0,1-((MAD*drSC)/(aSC*MADSC)).^2).^2;
        mcttd(MAD*drCT > aCT*MADCT,:) = [];
        mcctd(MAD*drCC > aCC*MADCC,:) = [];
        mst(MAD*drST > aST*MADST,:) = [];
        msc(MAD*drSC > aSC*MADSC,:) = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%   Recoard stats from this iteration   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    stats.Nclst(iter) = max(C);
    stats.Ne(iter)    = Ne;
    stats.Nct(iter)   = Nct;
    stats.Ncc(iter)   = Ncc;
    stats.Nst(iter)   = Nst;
    stats.Nsc(iter)   = Nsc;
    stats.Nsp(iter)   = Nsp;
    stats.DH(iter)    = s.DH;
    stats.DZ(iter)    = s.DH;
    stats.OS(iter,:)  = s.OS;
    stats.prctM(iter) = s.prctM;
    stats.iter(iter)  = s.iter;
    stats.jE{iter}    = jE;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%      Record final solution      %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hyp1(jE,2),hyp1(jE,1)] = m_xy2ll(1000*H(jE,1),1000*H(jE,2));
hyp1(jE,3) = H(jE,3);
jne = setdiff((1:stats.Ne0)',jE);
hyp1(jne,:) = NaN;

T1 = computeOriginTimesHypoTD(T0,H,S,phadir,jE,Ft,Fg,Fv,params.vpvs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      Bootstrap for error estimate          %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~params.Nboot || Ne < 2
    stats.err = [];
else
    
    % -- Final solution ray-tracing and
    % -- (need velocities and take-off angles) 
    [t,g,vto] = RayTrace1DhypoDD(H(uES(:,1),:),S(uES(:,2),:),Ft,Fg,Fv); 
    gCT  = [g(jptA1,:) g(jptB1,:) g(jptA2,:) g(jptB2,:)];
    gCC  = [g(jpcA1,:) g(jpcB1,:) g(jpcA2,:) g(jpcB2,:)];
    gST = [g(jpst1,:) g(jpst2,:)];
    gSC = [g(jpsc1,:) g(jpsc2,:)];
    vfCT = ones(Nct,1);
    vfCC = ones(Ncc,1);
    vfST = ones(Nst,1);
    vfSC = ones(Nsc,1);
    vfCT(jSct) = params.vpvs;
    vfCC(jScc) = params.vpvs; 
    vfST(jSst) = params.vpvs;
    vfSC(jSsc) = params.vpvs;
    vCT = [vto(jptA1)./vfCT vto(jptB1)./vfCT];
    vCC = [vto(jpcA1)./vfCC vto(jpcB1)./vfCC];
    vST = vto(jpst1)./vfST;
    vSC = vto(jpsc1)./vfSC;
    
    % -- Final differential times
    dtCT = mcttd(:,6);
    dtCC = mcctd(:,6);
    dtST = mst(:,5);
    dtSC = msc(:,5);
    
    % -- But using exact same equation set as final iteration    
    WCT = ones(Nct,1).*wdCT*params.wCT(end);
    WCC = ones(Ncc,1).*wdCC*params.wCC(end);
    WST = ones(Nst,1)*params.wST(end);
    WSC = ones(Nsc,1)*params.wSC(end);
    GCT = buildHypoTDmatrix(evsCT,gCT,vCT,WCT,Ne);
    GCC = buildHypoTDmatrix(evsCC,gCC,vCC,WCC,Ne);
    GST = buildStationDDmatrix(evsST,gST,vST,WST,Ne,1);
    GSC = buildStationDDmatrix(evsSC,gSC,vSC,WSC,Ne,1);
    dCT = WCT.*(dtCT-vfCT.*( (t(jptA1)-t(jptB1))-(t(jptA2)-t(jptB2)) ));
    dCC = WCC.*(dtCC-vfCC.*( (t(jpcA1)-t(jpcB1))-(t(jpcA2)-t(jpcB2)) ));
    dST = WST.*(dtST-vfST.*(t(jpst1)-t(jpst2)));
    dSC = WSC.*(dtSC-vfSC.*(t(jpsc1)-t(jpsc2)));
    
    % -- Note that no damping is included
    G = [GCT; GCC; GST; GSC; GZ];
    
    disp('--------------------------------------------------------------------')
    disp(['Bootstrapping (',num2str(params.Nboot),' iterations)'])
    stats.err = NaN*ones(stats.Ne0,3);
    stats.err(jE,:) = bootstrapHypoTDX(G,dCT,dCC,dST,dSC,dZ,H(jE,3),params);
end
disp(['hypoTDX complete ',datestr(datetime('now'))])
