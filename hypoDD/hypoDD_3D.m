function [hyp1,T1,stats] = hypoDD_3D(hyp0,T0,sta,phadir,mct,mcc,model,params)
% function [hyp1,T1,stats] = hypoDD_3D(hyp0,T0,sta,phadir,mct,mcc,sP,model,params)
%
% 2020-01-06
% This code mimics "hypoDD" by Felix Waldhauser.
% Some comments refer to Waldhauser and Ellsworth (2000,BSSA) as W2000
%
%   INPUTS 
%
%     hyp0 = [lat0,lon0,dep0] initial hypocenters (NEx3)
%       T0 = intial origin times (NEx1 datetime)
%      sta = [latS,lonS,depS(km)]
%   phadir = [EV,STA,PHA,T]
%      mct = catalog differential times [EVa,EVb,STA,PHA,DT,(CC)]
%            can be formed with "ph2dt.m"
%      mcc = waveform differential times [EVa,EVb,STA,PHA,DT,CC]
%            can be formed with "L1align" / "VCMC"
%    model = a structure with the important fields being:
%              V      = a 3D "interpolant" of velocity
%              mapcmd = a string of the command used to make the map
%                       projection for "V"
%              zrng   = minimum/maximum depths in km of the model.
%
%   params = structure of optional paramters
%
%
%  OUTPUTS
%
%     hyp1 = [lat,lon,dep] estimated hypocenter (NEx3)
%       T1 = estimated origin time (NEx1 datetime)

%    stats = a dictionary with msft,msft0...


if isempty(mct)
    mct = zeros(0,6);
end
if isempty(mcc)
    mcc = zeros(0,6);
end

Ne  = size(hyp0,1);
Ns  = size(sta,1);
Np  = size(phadir,1);
Nct = size(mct,1);
Ncc = size(mcc,1);

% -- Sort input/default parameters 
params = unpack_paramsHypoDD3D(params);

% -- EQ and Station coordinates in map units (km)
[ex,ey] = m_ll2xy(hyp0(:,2),hyp0(:,1));
[sx,sy] = m_ll2xy(sta(:,2),sta(:,1));
H = [ex/1000,ey/1000,hyp0(:,3)];
S = [sx/1000,sy/1000,sta(:,3)];

% -- MAD for gaussian noise (sigma_MAD in W2000 eq. 15)
MAD = 0.67449;

% -- Initialize model vector with origin times just set to zero
% -- Initialize vector of relocated events as full list
T1 = T0;
jE = (1:Ne)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    Intialize structure to record stats    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stats.Ne0 = Ne;
stats.Ns0 = Ns;
stats.Np0 = Np;
stats.Nct0 = Nct;
stats.Ncc0 = Ncc;

stats.Nclst = zeros(params.Niter,1);
stats.Ne    = zeros(params.Niter,1);
stats.Nct   = zeros(params.Niter,1);
stats.Ncc   = zeros(params.Niter,1);
stats.DH    = zeros(params.Niter,1);
stats.DZ    = zeros(params.Niter,1);
stats.DT    = zeros(params.Niter,1);
stats.OS    = zeros(params.Niter,4);
stats.prctM = zeros(params.Niter,1);
stats.jE    = cell(params.Niter,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    Print some initial stats    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['%%%%%%%%%%%%%%%   hypoDD starting ',datestr(datetime('now')),'   %%%%%%%%%%%%%%%'])
disp(['# earthquakes: ', num2str(Ne)])
disp(['# stations: ', num2str(Ns)])
disp(['# pick-based DDs: ',num2str(Nct)])
disp(['# waveform-based DDs: ',num2str(Ncc)])
disp(' ')
disp('Iter     Ne   Nclst      Nct      Ncc    DH(km)  DZ(km)  DT(s)  %d_msft')
disp('--------------------------------------------------------------------')

%fprintf('% 3d  % 5d   % 7d  % 7d  % 6.3f  % 6.3f  % 6.3f  % 5.2f\n',  ... 
%             iter,Ne,Nct,Ncc,sts.mnDH,sts.mnDV,sts.mnDT,sts.prctM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%      Begin DD iterations      %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:params.Niter
  
    % -- Cut-off thresholds of residuals
    aCT = params.alphCT(iter);
    aCC = params.alphCC(iter);
    dmax = params.dmax(iter); 
    lmbd = params.lmbd(iter);
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%      Cull events if necessary     %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -- Cluster events, ensuring continuous connection
    % -- Remove comparison for events in different clusters?
    % -- Unlike hypoDD, I do this separately for each iteration
    [mct,mcc,jE]   = checkMinStaHypoDD(mct,mcc,jE,params.minsta,params.minstaPS);
    [C,mct,mcc,jE] = clusterHypoDD(mct,mcc,jE,params.NobsC);
    
    Ne  = length(jE);
    Nct = size(mct,1);
    Ncc = size(mcc,1);
    
    % -- Check if we've eliminated all events!!
    if Ne < 2
        break
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%       Do Ray Tracing      %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % -- Assign "phase numbers" to each unique event/station combo
    uES = unique([mct(:,1) mct(:,3)
                  mct(:,2) mct(:,3)
                  mcc(:,1) mcc(:,3)
                  mcc(:,2) mcc(:,3)],'rows');
    uES(:,1) = uES(:,1);          
                  
    % -- Do Ray-Tracing
    [t,g,vto] = RayTrace3DhypoDD(H(uES(:,1),:),S(uES(:,2),:),model.V,... 
                                   params.minV,params.maxV);
 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Sort data prior to taking step %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,evsCT] = ismember(mct(:,1:2),jE);
    [~,evsCC] = ismember(mcc(:,1:2),jE);

    % -- Index is phase number
    [~,jpt1] =  ismember(mct(:,[1,3]),uES,'rows');  
    [~,jpt2] =  ismember(mct(:,[2,3]),uES,'rows');
    [~,jpc1] =  ismember(mcc(:,[1,3]),uES,'rows');
    [~,jpc2] =  ismember(mcc(:,[2,3]),uES,'rows');


    % -- Observed differential times
    dtCT = mct(:,5);
    dtCC = mcc(:,5);

    % -- Initialize weights from input (reweight P vs. S, CC vs. CT later)
    wCT  = mct(:,6);
    wCC  = mcc(:,6);

    % -- Indices of P, S equations, 
    jPct = find(mct(:,4)==1);
    jSct = find(mct(:,4)==2);
    jPcc = find(mcc(:,4)==1);
    jScc = find(mcc(:,4)==2);

    % -- Vector of Vp/Vs correction (for S)
    vfCT = ones(Nct,1);
    vfCC = ones(Ncc,1);
    vfCT(jSct) = params.vpvs;
    vfCC(jScc) = params.vpvs;

    % -- Distance between events, and corresponding weights
    wdCT = distanceWeightingHypoDD(mct,H,dmax);
    wdCC = distanceWeightingHypoDD(mcc,H,dmax);
    
    % -- Relative weighting of CT vs. CC  and P vs. S data
    WCT = wCT.*wdCT*params.wCT(iter);
    WCC = wCC.*wdCC*params.wCC(iter);
    WCT(jPct) = WCT(jPct)*params.wP(iter);
    WCT(jSct) = WCT(jSct)*params.wS(iter);
    WCC(jPcc) = WCC(jPcc)*params.wP(iter);
    WCC(jScc) = WCC(jScc)*params.wS(iter);

    % -- Index unique ray-paths to get take-off angles
    gCT = [g(jpt1,:) g(jpt2,:)];
    gCC = [g(jpc1,:) g(jpc2,:)];

    % -- Index unique ray-paths again to get velocities at source
    % -- Divide by Vp/Vs for S-data
    vCT = [vto(jpt1)./vfCT vto(jpt2)./vfCT];
    vCC = [vto(jpc1)./vfCC vto(jpc2)./vfCC];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  Actually build linear system and solve %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % -- Form core DD matrices and RHS vectors 
    GCT = buildHypoDDmatrix(evsCT,gCT,vCT,WCT,Ne);
    GCC = buildHypoDDmatrix(evsCC,gCC,vCC,WCC,Ne);
    dCT = WCT.*(dtCT-vfCT.*(t(jpt1)-t(jpt2)));
    dCC = WCC.*(dtCC-vfCC.*(t(jpc1)-t(jpc2)));
    
    % -- Form zero-shift  and damping matrices and RHS vectors
    [GZ,dZ] = ZeroShiftHypoDD(C,params.W00);
    [GD,dD] = DampingMatrixHypoDD(C,stats.Ne0,lmbd);
    
    % -- Combine portions of systems and solve with Newton's method
    G = [GCT; GCC; GZ; GD];
    d = [dCT; dCC; dZ; dD];
    
    [m,s] = solveHypoDD(G,d,H(jE,3),params.NewtonSteps,params.CGsteps, ... 
                                    params.minZ,params.maxZ);
                              
    fprintf(' % 3d   % 5d   % 3d   % 7d  % 7d   % 6.3f  % 6.3f  % 6.3f  % 5.2f\n',  ... 
             iter,Ne,max(C),Nct,Ncc,s.DH,s.DZ,s.DT,s.prctM)                           

    % -- Update hypocenter array with computed optimal shift
    H(jE,1) = H(jE,1) + m(1:Ne);
    H(jE,2) = H(jE,2) + m(Ne+[1:Ne]);
    H(jE,3) = H(jE,3) + m(2*Ne+[1:Ne]);
    T1(jE)  = T1(jE)  + seconds(m(3*Ne+(1:Ne)));
    
    % -- Update differential times...
    mct(:,5) = mct(:,5) + m(3*Ne+evsCT(:,2)) - m(3*Ne+evsCT(:,1));
    mcc(:,5) = mcc(:,5) + m(3*Ne+evsCC(:,2)) - m(3*Ne+evsCC(:,1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Check for outlier equations %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iter < params.Niter
        % -- Indices, Weights for CT,CC,...
        jct = (1:Nct)';
        jcc = Nct+(1:Ncc)';

        % -- Reweight and cull equations based on residual (eq. 15 of W2000)
        drCT  = G(jct,:)*m-d(jct);
        drCC  = G(jcc,:)*m-d(jcc);
        MADCT = median(abs(drCT-median(drCT)));
        MADCC = median(abs(drCC-median(drCC)));

        mct(:,6) = max(0,1-((MAD*drCT)/(aCT*MADCT)).^2).^2;
        mcc(:,6) = max(0,1-((MAD*drCC)/(aCC*MADCC)).^2).^2;
        mct(MAD*drCT > aCT*MADCT,:) = [];
        mcc(MAD*drCC > aCC*MADCC,:) = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%   Recoard stats from this iteration   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    stats.Nclst(iter) = max(C);
    stats.Ne(iter)    = Ne;
    stats.Nct(iter)   = Nct;
    stats.Ncc(iter)   = Ncc;
    stats.DH(iter)    = s.DH;
    stats.DZ(iter)    = s.DH;
    stats.DT(iter)    = s.DT;
    stats.OS(iter,:)  = s.OS;
    stats.prctM(iter) = s.prctM;
    stats.jE{iter}    = jE;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%      Record final solution      %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hyp1(jE,2),hyp1(jE,1)] = m_xy2ll(1000*H(jE,1),1000*H(jE,2));
hyp1(jE,3) = H(jE,3);

jne = setdiff((1:stats.Ne0)',jE);
hyp1(jne,:) = NaN;
T1(jne) = NaT;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      Bootstrap for error estimate          %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~params.Nboot
    stats.err = [];
else
    
    % -- Final solution ray-tracing and
    % -- (need velocities and take-off angles) 
    [t,g,vto] = RayTrace3DhypoDD(H(uES(:,1),:),S(uES(:,2),:),model.V,... 
                                   params.minV,params.maxV);       
    gCT  = [g(jpt1,:) g(jpt2,:)];
    gCC  = [g(jpc1,:) g(jpc2,:)];
    vfCT = ones(Nct,1);
    vfCC = ones(Ncc,1);
    vfCT(jSct) = params.vpvs;
    vfCC(jScc) = params.vpvs;
    vCT = [vto(jpt1)./vfCT vto(jpt2)./vfCT];
    vCC = [vto(jpc1)./vfCC vto(jpc2)./vfCC];
    
    % -- Final differential times
    dtCT = mct(:,5);
    dtCC = mcc(:,5);
    
    % -- But using exact same equation set as final iteration                           
    WCT = ones(Nct,1)*params.wCT(end);
    WCC = ones(Ncc,1)*params.wCC(end);
    GCT = buildHypoDDmatrix(evsCT,gCT,vCT,WCT,Ne);
    GCC = buildHypoDDmatrix(evsCC,gCC,vCC,WCC,Ne);
    dCT = WCT.*(dtCT-vfCT.*(t(jpt1)-t(jpt2)));
    dCC = WCC.*(dtCC-vfCC.*(t(jpc1)-t(jpc2)));
    
    G = [GCT; GCC; GZ];
    
    disp('--------------------------------------------------------------------')
    disp(['Bootstrapping (',num2str(params.Nboot),' iterations)'])
    stats.err = NaN*ones(stats.Ne0,4);
    stats.err(jE,:) = bootstrapHypoDD(G,dCT,dCC,dZ,H(jE,3),params);
end
disp(['hypoDD complete ',datestr(datetime('now'))])
