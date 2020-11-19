function [hyp,T,stats] = hypoTD(hyp0,T0,stadir,phadir,mcc,sP,model,params)
% function [hyp,T,stats] = hypoTD(hyp0,T0,stadir,phadir,mcc,sP,model,params)
%
% 2020-05-01
% This function is my own "triple-difference" relative hypocenter
% inversion code. The triple-difference is:
%
%  ((T1a-T2a)-(T1b-T2b))_obs - ((T1a-T2a)-(T1b-T2b))_cal
%   (stations 1,2, events a,b0
%     ...need to compute derivatives?? ...or use grid search somehow??
%
% by using the TD instead of DD, we ignore the origin time throughtout the
% main inversion---solving for a best-fit origin time at the end.
%
%  % Each step
%  %  1) compute ray take-off directions for each phase
%  %  2) compute derivates dTia/dx
%
%   The solution vector "x" has
%   The matrix "A" has all  
%
%
%   INPUTS 
%
%     hyp0 = [lat0,lon0,dep0] initial hypocenters (NEx3)
%       T0 = intial origin times (NEx1 datetime)
%   stadir = [latS,lonS]
%   phadir = [EV,STA,PHA,T]
%      mcc = modified "mcctimes" matrix [EVa,EVb,STA,PHA,DT,CC]
%       sP = [EV,STA,TsP]
%    model = model name 'GSC','ALASKA','ma2011'
%   params = structure of optional paramters
%
%  OUTPUTS
%
%      hyp = [lat,lon,dep] estimated hypocenter (NEx3)
%        T = estimated origin time (NEx1 datetime)

%    stats = a dictionary with msft,msft0...
%%%%%%%%%%%%%%%%%%%%%%%
%    picks = [[picks], res (Tpred-Tobs), D epicentral, azimuth, take-off angle]
%    picks = {[latS,lonS,1/2 for P/S,T (relative to T0 in s)]}
%             it may also have columns for [res,dist_epi,azi,toa]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- SET UP
unpack_paramsHypoTD;

Ne  = size(hyp0,1);
Ns  = size(stadir,1);
Np  = size(phadir,1);
NsP = size(sP,1);


% -- Prepare grids of travel-times, take-off angles, velocity
prep_grid_hypoTD;

% ((T1a-T2a)-(T1b-T2b))_obs - ((T1a-T2a)-(T1b-T2b))_cal
%  == ( (T1a-T1b)-(T2a-T2b) )_obs  - ( (T1a-T1b)-(T2a-T2b) )_cal
% -- Need to find triple differences from mcc
% -- Could also do equivalent for dt.ct comparisons...
% -- First should reduce mcc to [EVa,EVb,STA,PHA,DT,CC,Pa,Pb] 
% -- (average over instruments)
% -- Also must create "dt.ct","catalog" version of mcc
%      mct = "dt.ct" verstion [EVa,EVb,STA,PHA,DT,CC]

evC = SelectEventCombinations(hyp0,phadir,max(max(mxDCT),max(mxDCC)),mxNB);
mct = phadir2dt(phadir,hyp0,evC);

% -- Epicentral distance, azimuth, from (Event x Station)
[dES,aES] = distance(repmat(hyp0(:,1),1,Ns),repmat(hyp0(:,2),1,Ns), ...
                     repmat(stadir(:,1)',Ne,1),repmat(stadir(:,2)',Ne,1),[6371,0]);
zES = repmat(hyp0(:,3),1,Ns);     

% -- Predict times, angles, by interpolation of previous table 
thtES = interp2(Zrng,Drng,tht,zES,dES);


jct = sub2ind([Ne Ns],mct(:,1),mct(:,3));
jcc = sub2ind([Ne Ns],mcc(:,1),mcc(:,3));

gESct = [cosd(aES(jct)).*sind(thtES(jct)) ...
         cosd(aES(jct)).*cosd(thtES(jct)) ...
         cosd(thtES(jct))];
gEScc = [cosd(aES(jcc)).*sind(thtES(jcc)) ...
         cosd(aES(jcc)).*cosd(thtES(jcc)) ...
         cosd(thtES(jcc))];
     
mttcc = TripleDifference(mcc,gEScc,mxTD); 
mttct = TripleDifference(mct,gESct,mxTD);

stats.NTDcc = size(mttcc,1);
stats.NTDct = size(mttct,1);
disp(['Triple-Differences Complete  ',datestr(datetime('now'))])

hyp = hyp0;
hyp(hyp(:,3) < minZ,3) = minZ;
hyp(hyp(:,3) > maxZ,3) = maxZ;

% -- Convert lat/lon to map coorditnates (in km)
[ex,ey] = m_ll2xy(hyp(:,2),hyp(:,1));
mx = mean(ex);
my = mean(ey);
ex = (ex-mx)/1e3;
ey = (ey-my)/1e3;
ez = hyp(:,3);
x00 = [ex; ey; ez];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0;
msft = zeros(max_iter+1,6);

NcullCC = 0;
NcullCT = 0;
NcullSP = 0;

for iter = 1:max_iter
    
    % -- Prepare for iteration
    prep_iteration_hypoTD;
    build_ACC;
    build_ACT;    
    build_AsP;
    build_A0;
    
    A = [wCC*WCC*ACC; wCT*WCT*ACT; wsP*AsP; w0*A0];
    b = [wCC*WCC*bCC; wCT*WCT*bCT; wsP*bsP; w0*b0];
    
    if iter==1
        msft(1,:) = [norm(A*x0-b), norm(WCC*(ACC*x0-bCC)), norm(WCT*(ACT*x0-bCT)),...
                     norm(AsP*x0-bsP), norm(A0*x0-b0)]; 
    end
    
    % -- Solve! I should explore different options here
    
    % -- It would be nice to use L1...
    % -- If I used L2 + steepest descent method I could use actual direct
    % -- picks instead of x0, x00 matrices
    x = solveHypoTD(A,b,x0,20,20,minZ,maxZ);
    disp(['Iter ',num2str(iter),'/',num2str(max_iter),'  mean(dx) = ', ...
          num2str(mean(abs(x-x0))),'  ',datestr(datetime('now'))])
    
    % -- Convert x,y back to lat/lon, update hypocenters
    ex = 1e3*x(1:Ne)+mx;
    ey = 1e3*x(Ne+[1:Ne])+my;
    
    [hyp(:,2),hyp(:,1)] = m_xy2ll(ex,ey);
    hyp(:,3) = x(2*Ne+[1:Ne]);
      
    % -- Compute misfit
    msft(iter+1,:) = [norm(A*x-b), norm(WCC*(ACC*x-bCC)), norm(WCT*(ACT*x-bCT)),...
                      norm(AsP*x-bsP), norm(A0*x-b0)]; 
    
    
    % -- ADD OUTLIER REMOVAL!!!   
    RemoveOutliersHypoTD;
    
    if mean(abs(x-x0)) < 1e-5
        break
    end
end


% -- Compute new origin times based on final hypocenter
prep_iteration_hypoTD;
computeOriginTimes;

% -- Relative weight of each part of matrix (CT,CC,sP,A0,A00)
% -- For each event... sum norms of equation coefficients?
RecordMatrixWeighting;


% -- Organize output stats
stats.msft     = msft;
stats.mtrxW    = mtrxW;
stats.relMtrxW = relMtrxW;
stats.wghts    = [wCC wCT wsP w0 w00 wMZ];
stats.NcullCC  = NcullCC;
stats.NcullCT  = NcullCT;
stats.NcullSP  = NcullSP;
stats.x        = x;
stats.params   = params;
stats.iter     = iter;

disp(['solveHypoTD Complete  ',datestr(datetime('now'))])
