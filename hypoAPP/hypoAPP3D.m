function [H,T,P,S] = hypoAPP3D(H,T,P,V,params)
% function [H,T,P,S] = hypoAPP3D(H,T,P,V,params)
%
% 2020-03-14
% This function estimates an earthquake hypocenter by in a 3D velocity
% model using an oct-tree iterative grid-search method from NonLinLoc, 
% which is described concisely here: 
%
%        http://alomax.free.fr/nlloc/octtree/OctTree.html
% 
% or in Lomax et al. (2000) "Probabilistic Earthquake Location in 3D
% and Layered Models"
%
%
%   INPUTS 
%
%        H = initial hypocenter estimate [x,y,z] (km) z==depth
%        T = intial origin time estimate(datetimes)
%        P = Pick matrix: [xS,yS,zS,PS W,tObs]
%            zS   == Station depth!! Note not eleveation!!(km)
%            PS   == 1/2 for P/S pick
%            W    == relative weight of pick
%            tObs == observed travel time (s)
%        V = griddedInterpolant or scatteredInterpolant of velocity
%            with x/y/z coordinates in km, in the same map projection 
%            used to define the H coordinates. I use the m_map package
%            to define a cartersian system, but divide horizonatal 
%            coordinates by 1000 to convert from m to km. I also generally
%            apply an earth-flattening transformation to the vertical 
%            coordinates and velocities befor makeing the interpolant.
%
%   params = an (optional) structure with the following fields:
%         vpvs = constant Vp/Vs applied to whole model, default is sqrt(3)
%          tol = minimum step length to worry about (km)
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
%         minV = minimum bound on P-velocity
%         maxV = maximum bound on P-velocity
%         sclF = If the cube being divided has side length <= sclF (km)
%                hypoAPP3D will only call RayTrace3D once, and just 
%                perturb the ray-path 
%         
%
%  OUTPUTS
%
%       H = Final hypocenter [x,y,z] (km) z==depth 
%       T = Final origin time (datetimes)
%       P = Pick matrix: [xS,yS,zS,PS W,tObs,tPred,res]
%       S = a "structure" with msft,msft0...


if nargin < 5
    params = struct();
end

% -- Sort input/default parameters (save initial step length)
params   = unpack_paramsHypoAPP3D(params);

% -- Check that depth is within bounds
H( H(:,3) < params.minZ ,3) = params.minZ;
H( H(:,3) > params.maxZ ,3) = params.maxZ;

% -- Default "exit" code is for max. iterations reachs
S.exit = "Max_Iter";

% -- Copy params over to stats stored in "S"
field = fieldnames(params);
for ii = 1:length(field)
    S.(field{ii}) = params.(field{ii});
end
  
% -- Initialize other stats fields 
S.NP    = 0;
S.NS    = 0;
S.dof   = 0;
S.iter  = 0;
S.res0  = [];
S.res   = [];
S.msft0 = Inf;
S.msft  = Inf;
S.errX  = Inf*ones(1,2);
S.errY  = Inf*ones(1,2);
S.errZ  = Inf*ones(1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Prep event
Np       = size(P,1);

% -- Get unique stations...ray-tracing is the same for P/S to same station
x1  = [P(:,1:2), P(:,3)];
ux1 = unique(x1,'rows');
[~,jsta] = ismember(x1,ux1,'rows');
Nus = size(ux1,1);

% -- Record initial misfit, check for stations with ray-tracing issues...
tpp   = RayTrace3D(repmat(H,Nus,1),ux1,V,params.minV,params.maxV);
PSfac = (params.vpvs).^(P(:,4)-1);
tPred = tpp(jsta).*PSfac;
jnan  = find(isnan(tpp(jsta)));

if length(jnan)
    P(jnan,:)   = [];
    PSfac(jnan) = [];
    tPred(jnan) = [];
    Np = size(P,1);
    
    % -- Check if there are still a sufficient number of stations
    % -- Note that this code assumes there were enough to begin with
    if Np < params.minpha || ~prod(ismember([1;2],P(:,4)))
        S.exit = 'MinPha_NaN';
        return
    end   
end

% -- Normalize pick weights and label observed pick times
% -- Compute effective degrees of freedom (dof)
wght  = Np*P(:,5)/sum(P(:,5));
dof   = Np^2/sum(wght.^2);
tObs  = P(:,6); 


% -- Check for outlier residuals?
S.res0 = (tObs-tPred)-WeightedMedian(tObs-tPred,wght);
[mxres,jmx] = max(abs(S.res0));
while mxres > params.picktol*median(abs(S.res0))
    P(jmx,:)    = [];
    tPred(jmx)  = [];
    PSfac(jmx)  = [];
    tObs(jmx)   = [];
    Np          = Np-1;
    wght        = Np*P(:,5)/sum(P(:,5));
    S.res0      = (tObs-tPred)-WeightedMedian(tObs-tPred,wght);
    [mxres,jmx] = max(abs(S.res0));
end

% -- Check if there are still a sufficient number of stations
% -- Note that this code assumes there were enough to begin with
if Np < params.minpha || ~prod(ismember([1;2],P(:,4)))
    S.exit = 'MinPha_Tol';
    return
end
    
% -- Reset station list if any were culled
x1  = [P(:,1:2), P(:,3)];
ux1 = unique(x1,'rows');
[~,jsta] = ismember(x1,ux1,'rows');
Nus = size(ux1,1);

% -- Compute initial misfit
% -- Define some scale to be applied to the misfits before converting
% -- to probabiliities in order to keep proababilites as reasonable orders
% -- of magnitude
S.msft0 = sum(wght.*abs(S.res0));
Mscl    = S.msft0.^-1;


% -- Iterations will divide one cell into eight--prepare grid for this
[xg,yg,zg] = ndgrid(-0.5:0.5,-0.5:0.5,-0.5:0.5);
xg = xg(:);
yg = yg(:);
zg = zg(:);
NnD = length(xg);


% -- Create grid points. Ensure search volume respects minZ/maxZ
% -- What to do when Lcube/2 > depth H(3)?
scl0 = params.Lcube/3;
zgC  = min( max(H(3),params.minZ+params.Lcube/2), ... 
                     params.maxZ-params.Lcube/2);
[xg0,yg0,zg0] = ndgrid(-1:1,-1:1,-1:1);
xg0 = scl0*xg0(:)+H(1);
yg0 = scl0*yg0(:)+H(2);
zg0 = scl0*zg0(:)+zgC;
Ng0 = length(xg0);

% -- Store irregular grid [xg,yg,zg,scl,prob]
% --  xg,yg,zg define grid point centers (of cubic volume)
% --  scl is the side length of the cube
% --  prob is the unscaled probaility of that cube (msft.^(-Np))
A  = [xg0 yg0 zg0 scl0*ones(Ng0,1) zeros(Ng0,1)];

% -- More points should be sampled in shallow region if
% -- intial hypocenter was shallower than Lcube/2...
% -- (Divide them into 8?) This increases the number of initial points
% -- from 27 to 90...
if zgC > H(3)
    scl = scl0/2;
    zm = min(A(:,3));
    for jj = find(A(:,3)==zm)'
        An = [A(jj,1)+scl*xg, A(jj,2)+scl*yg, A(jj,3)+scl*zg, scl*ones(NnD,1) zeros(NnD,1)];
        A = [A; An];
    end
    A(A(:,3)==zm,:) = [];
end
Ng  = size(A,1);

for ii = 1:Ng
    x0      = repmat(A(ii,1:3),Nus,1);
    tpp     = RayTrace3D(x0,ux1,V,params.minV,params.maxV);
    tPred   = tpp(jsta).*PSfac;
    res     = (tObs-tPred)-WeightedMedian(tObs-tPred,wght);
    msft    = sum(wght.*abs(res));
    A(ii,5) = (A(ii,4)^3)*(msft*Mscl).^(-dof);
end

% -- Create matrix to store previous hypocenters
% -- Compute hypocenter as centroid of probability
% -- For efficiency compute weighted mean of probability, rather
% -- than weighted median
H      = [H; H; Inf*ones(params.NstepTol-1,3)];
w      = A(:,5)/sum(A(:,5));
H(1,:) = sum([w.*A(:,1), w.*A(:,2), w.*A(:,3)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Rest of iterations
for iter = 1:params.max_iter

    % -- Initialize previous hypocenters
    H(2:end,:) = H(1:end-1,:);
    
    % -- Choose cell to divide into eight
    [~,jmx] = max(A(:,5));
    scl = A(jmx,4)/2;
    
    % -- Divide cell into eight and do ray-tracing, compute
    % -- misfits/probabilities for each
    An  = [A(jmx,1)+scl*xg, A(jmx,2)+scl*yg, A(jmx,3)+scl*zg, scl*ones(NnD,1) zeros(NnD,1)];
    Nn  = size(An,1);
    vol = scl^3;
    
    % -- If the cube being divided is large, trace rays separately for 
    % -- each new point.
    if scl > params.sclF  
        for ii = 1:Nn
            x0       = repmat(An(ii,1:3),Nus,1);
            tpp      = RayTrace3D(x0,ux1,V,params.minV,params.maxV);
            tPred    = tpp(jsta).*PSfac;
            res      = (tObs-tPred)-WeightedMedian(tObs-tPred,wght);
            msft     = sum(wght.*abs(res));
            An(ii,5) = vol*(msft*Mscl).^(-Np);
        end 
        
    % -- If the cube is small, don't do full ray-tracing.
    % -- Assume ray takes basically the same shape, just resample the 
    % -- velocity model with the slightly pertubed ray-path
    else
        x0     = repmat(A(jmx,1:3),Nus,1);
        [~,RP] = RayTrace3D(x0,ux1,V,params.minV,params.maxV);
        wdxrp  = cell(Nus,1);
        for iss = 1:Nus
            wdxrp{iss} = linspace(1,0,length(RP{iss}))';
        end
        tpp = zeros(Nus,1);
        for ii = 1:Nn      
            for iss = 1:Nus
                rp = RP{iss}+wdxrp{iss}*(An(ii,1:3)-A(jmx,1:3));
                Vt = Vfun(rp,V,params.minV,params.maxV);
                di = sqrt(sum(diff(rp).^2,2));
                tpp(iss) = 0.5*sum( di.*(1./Vt(1:end-1) + 1./Vt(2:end)) );
            end
            tPred    = tpp(jsta).*PSfac;
            res      = (tObs-tPred)-WeightedMedian(tObs-tPred,wght);
            msft     = sum(wght.*abs(res));
            An(ii,5) = vol*(msft*Mscl).^(-Np);
        end
    end
     
    % -- Merge new cell with existing ones
    % -- Remove the cell that was just divided into eight
    A = [A([1:jmx-1,jmx+1:end],:); An];
    
    % -- Compute new hypocenter and check if converged
    w = A(:,5)/sum(A(:,5));
    H(1,:) = sum([w.*A(:,1), w.*A(:,2), w.*A(:,3)]);
    
    dH = norm(H(1,:)-H(end,:));
    if norm(H(1,:)-H(end,:)) < params.tol
        S.exit = 'tol';
        break
    end
    
end % -- Iterations Loop


% -- Redo last step to get uncertainy bounds
H = H(1,:);
        
% -- Trace rays for final hypocenter estimate
x0     = repmat(H,Nus,1);
tpp    = RayTrace3D(x0,ux1,V,params.minV,params.maxV);
tPred  = tpp(jsta).*PSfac;

% -- Compute updated origin time, adjust tObs according, and compute
% -- final residuals and misfit
dt0    = WeightedMedian(tObs-tPred,wght);
tObs   = tObs-dt0;
T      = T+seconds(dt0);
S.res  = tObs-tPred;
S.msft = sum(wght.*abs(S.res));

% -- Build final pick matrix [xS,yS,zS,PS W,tObs,tPred,res]
P = [P(:,1:4), wght, tObs, tPred, res];

% -- Redo last step to get uncertainy bounds
% -- Uncertainties (95%?) in x/y/z directions
pct = [(1-params.CI)/2  1-(1-params.CI)/2];
[~,S.errX] = WeightedMedian(A(:,1),A(:,5),pct);
[~,S.errY] = WeightedMedian(A(:,2),A(:,5),pct);
[~,S.errZ] = WeightedMedian(A(:,3),A(:,5),pct);
S.errX = S.errX-H(1);
S.errY = S.errY-H(2);
S.errZ = S.errZ-H(3);

% -- Record some other stats / parameters used
S.iter     = iter;
S.NP       = length(find(P(:,4)==1));
S.NS       = length(find(P(:,4)==2));
S.dof      = dof;
    


% -- Divide center point right away? (jC = 14th of 27 initial points)
%jC  = fix((Ng0-1)/2)+1;
%A = [A; A(jC,1)+scl*xg, A(jC,2)+scl*yg, A(jC,3)+scl*zg, scl*ones(NnD,1) zeros(NnD,1)];
%A(jC,:) = [];

%jj = [1:jmx-1,jmx+1:size(A,1)];
%figure(1)
%plot3(A(:,1),A(:,2),A(:,3),'bo')
%plot3(A(jj,1),A(jj,2),A(jj,3),'bo',An(:,1),An(:,2),An(:,3),'ro',A(jmx,1),A(jmx,2),A(jmx,3),'go')
%axis equal
%keyboard
    
%H(1,:) = [WeightedMedian(A(:,1),A(:,5)), ...
%          WeightedMedian(A(:,2),A(:,5)), ... 
%          WeightedMedian(A(:,3),A(:,5))];

%An(An(:,3) < params.minZ,:) = []; This is unnecessary.
%An(An(:,3) > params.maxZ,:) = [];
