function [TRND,PLNG,R,uS,uR] = StressInvV2014(sdr,nx,ndg,MTHD,prcntU)
% function [TRND,PLNG,R,uS,uR] = StressInvV2014(sdr,nx,ndg,MTHD,prctU)
%
% 2019-11-27
% This function mainly follows the Vavrycuk [2014,GJI] iterative method 
% for regional stress inversion for focal mechanisms with unknown true
% fault planes. Each of the 'nx' iteration perturbs the focal mechanisms 
% random amounts according to 'ndg', and then tries to solve for a stress
% tensor while choosing the correct fault planes (an iterative process in
% itself). Several methods of choosing the 'correct' fault plane are 
% available via the input 'MTHD'
%
% 2020-07-13
% Updated to perturb MTs by rotating the whole focal sphere, instead of
% altering strike/dip/rake. Also it now allows for errors speciifc to 
% each event
%
% 2021-12-21
% Updated to rotate focal mechanisms (to perturb them according to their
% uncertainty) by using a random Quaternion
%
%    INPUTS
%
%       sdr = NE x 6 matrix [STRIKE1,DIP1,RAKE1,STRIKE2,DIP2,RAKE2]
%        nx = number of iterations (random perturbations)
%       ndg = assumed error on focal mechanisms (in degrees)
%             can be NE x 1 vector. e.g. a standard deviation
%             or equivalent L1 measure
%      MTHD = 'method' can be:
%             'I' = optimize fault instability
%             'R' = optimize predicted vs. observed rakes
%             'C' = optimize a combination of fault instability and 
%                    rake according to I.*( 1 + dot(r_obs,r_pre) )
%                    'C' is the default
%    prcntU = percentiles to return as uncertainty bounds
%             can be e.g. [68,95] (which is the default)
%
%   OUTPUTS
%
%     TRND = "trends" of principal stress directions [sg1,sg2,sg3]
%            (horizontal angle in degrees from north)
%     PLNG = "plunges" of principal stress directions [sg1,sg2,sg3]
%        R = stress ratio (sg1-sg2)/(sg1-sg3). 
%       uS = uncertainty of principal stress directions in degrees,
%            corresponds to percentiles given in prctU
%       uR = uncertainty in stress ratio corresponding to prcntU
%            computed as prctile(abs(Rs-Rfinal),prcntU)

if nargin < 5
    prcntU = [68 95];
    if nargin < 4
        MTHD = 'C';
    end
end

% -- Fault friction (mu), one for each run of experiment
% -- Noise levels (ndg), standard deviation error on MTs (in degrees)
mus = 0.45+0.44*(2*(rand(nx,1)-0.5));
ne  = size(sdr,1);
if length(ndg)==1
    ndg = ndg*ones(ne,1);
end


% -- decide what function to maximize
if strcmp(MTHD,'I')
    Qfun = @(I,dotR) I;
elseif strcmp(MTHD,'R')
    %Qfun = @(I,dr) -abs(dr);
    Qfun = @(I,dotR) dotR;
elseif strcmp(MTHD,'C')
    %Qfun = @(I,dotR) I.*(1-abs(dotR)/90);
    Qfun = @(I,dotR) I.*(dotR);
end

tsdr = zeros(ne,3);
tn = zeros(ne,3);
tv = zeros(ne,3);

% -- Array to store output stress tensors (in vector form), msft
xf   = zeros(5,nx);
msft = zeros(nx,1);
cnvg = zeros(nx,1);
fpf  = zeros(ne,nx);
Vs   = zeros(3,3,nx);
Rs   = zeros(nx,1);


[n0,v0] = sdr2nv(sdr);
n1 = n0;
v1 = v0;

% -- Run in parallel if pool is open
pool = gcp('nocreate');
if isempty(pool)
    Nworker = 0;
else
    Nworker = pool.NumWorkers;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for ii = 1:nx
parfor (ii=1:nx,Nworker)

    mu = mus(ii);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Add noise to FMs
    
    % -- Angle to rotate each event (in degrees)
    %NdgR = ndg.*randn(ne,1);
    
    % -- Random Laplacian (with "ndg" == b, ... variance=2*b^2);
    u = rand(ne,1)-0.5;
    NdgR = abs(ndg.*sign(u).*log(1-2*abs(u)));
    
    for ie = 1:ne
        
        % -- Random Quaternion to perform rotation (perturbations)
        RM = Quat2RotMtrx(RandomQuaternion(1));
        RM = real(RM^(NdgR(ie)/acosd((trace(RM)-1)/2)));
        
        n1(ie,:) = n0(ie,:)*RM;
        v1(ie,:) = v0(ie,:)*RM;
    end
    
    % -- Not sure why the results are coming out with a (zero) imaginary
    % -- component...
    n1 = real(n1);
    v1 = real(v1);
    n2 = v1;
    v2 = n1;
        
    A1 = BuildStressMatrix(n1);
    A2 = BuildStressMatrix(n2);
    
    rn = zeros(ne,3,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % -- Start by randomly selecting fault planes (assign 1 or 2)
    fp = rand(ne,1); 
    fp(fp>0.5) = 2;
    fp(fp<2)   = 1;
           
    for iter = 1:5
              
        tn = n1;
        tv = v1;
        jj2 = find(fp==2);
        tn(jj2,:) = n2(jj2,:);
        tv(jj2,:) = v2(jj2,:);
        
        A = BuildStressMatrix(tn);
        b = -tv(:);
                
        x = A\b;
        %[x,msft(ii)] = iterativeStressSolver(A,b,1e-4);
        
        X = x([1 2 3; 2 4 5; 3 5 1]);
        X(3,3) = -X(1,1)-X(2,2);

        % -- Decompose the stress tensor, ensure eigenvalues go from most 
        % -- positive (compressional) to most negative (tensional)
        % -- Also ensure eigenvetors point down.
        [V,S] = eig(X);
        S = flip(diag(S));
        V = fliplr(V);
        R = (S(1)-S(2))/(S(1)-S(3));  

    
        % -- Need to rotate normal to be relative to stress axes, then compute
        % -- instability (via sigma,tau,mu)
        % -- This follows Vavrycuk [2014,GJI] eqs. 16 to 18
        rn(:,:,1) = n1*V;
        rn(:,:,2) = n2*V;
        sig = squeeze(rn(:,1,:).^2 + (1-2*R)*rn(:,2,:).^2 - rn(:,3,:).^2);
        tau = sqrt(squeeze(rn(:,1,:).^2 + (1-2*R)^2*rn(:,2,:).^2 + rn(:,3,:).^2)-sig.^2); 
        I   = (tau - mu*(sig-1)) / (mu + sqrt(1+mu^2));
        
        
        % -- Predict slip vectors, normalize them to length 1       
        pv1 = reshape(-A1*x,[ne,3]);
        pv2 = reshape(-A2*x,[ne,3]);
        pv1 = pv1./repmat(vecnorm(pv1')',1,3);
        pv2 = pv2./repmat(vecnorm(pv2')',1,3);
        
        % -- Take dot product of predicted and actual slip vectors
        dotR = [sum(pv1.*v1,2), sum(pv2.*v2,2)];
        
        Q = Qfun(I,dotR);
        [~,fpn] = max(Q,[],2);       
        
        
        if sum(fp==fpn)==ne 
            fp = fpn;  
            cnvg(ii) = 1;
            break
        else
            fp = fpn;
        end    
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xf(:,ii)   = x;
    %msft(ii)   = norm(A*x-b,1)/norm(b,1);
    fpf(:,ii)  = fpn;
    Rs(ii)     = R;
    Vs(:,:,ii) = V;
    
end


% -- Take a weighted mean of all output stress tensors
jc  = find(cnvg);
njc = length(jc);

% -- Get axes estimates from mean X...
xmn = sum(xf,2)/nx;
X      = xmn([1 2 3; 2 4 5; 3 5 1]);
X(3,3) = -X(1,1)-X(2,2);

[TRND,PLNG,R] = StressTensor2Axes(X);

V = zeros(3);
V(:,1) = [cosd(90-TRND(1))*cosd(PLNG(1)); sind(90-TRND(1))*cosd(PLNG(1)); sind(PLNG(1))];
V(:,2) = [cosd(90-TRND(2))*cosd(PLNG(2)); sind(90-TRND(2))*cosd(PLNG(2)); sind(PLNG(2))];
V(:,3) = [cosd(90-TRND(3))*cosd(PLNG(3)); sind(90-TRND(3))*cosd(PLNG(3)); sind(PLNG(3))];

% -- Final sigma 1,2,3 direction are in the columns of V
% -- each iteration's version are in Vs..
NU = length(prcntU);
pU = prcntU(:);

xyz1 = squeeze(Vs(:,1,:));
xyz2 = squeeze(Vs(:,2,:));
xyz3 = squeeze(Vs(:,3,:));
da1 = acosd(V(:,1)'*xyz1)';
da2 = acosd(V(:,2)'*xyz2)';
da3 = acosd(V(:,3)'*xyz3)';
da1(da1>90) = 180-da1(da1>90);
da2(da2>90) = 180-da2(da2>90);
da3(da3>90) = 180-da3(da3>90);

sg1 = prctile(da1,pU);
sg2 = prctile(da2,pU);
sg3 = prctile(da3,pU);
uS  = [sg1,sg2,sg3];

bU = 50+[-pU'; pU']/2;
uR = reshape(prctile(Rs,bU(:)),[2 NU]);



