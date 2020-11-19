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
% 2020-09-24
% This was copied from StressInvV2014, then I added another part where
% it makes contour plots of the three principal stress axes probabilities.
%
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
    Qfun = @(I,dotR) dotR;
elseif strcmp(MTHD,'C')
    Qfun = @(I,dotR) I.*(dotR);
end


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:nx

    mu = mus(ii);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Add noise to strk,dip,rake 
    
    % -- Angle to rotate each event (in degrees)
    %NdgR = ndg.*randn(ne,1);
    
    % -- Random Laplacian (with "ndg" == b, ... variance=2*b^2);
    u = rand(ne,1)-0.5;
    NdgR = -ndg.*sign(u).*log(1-2*abs(u));
    
    for ie = 1:ne
        
        % -- Pythagorean theorem applies to rotations around perperdicular
        % -- axes?... sepearate rotation of NdgR deg into x,y,z components
        % -- IS THIS TRUE???
        a = randn(3,1);
        a = a*NdgR(ie)/norm(a);
        RMx = [1 0 0; 0 cosd(a(1)) -sind(a(1)); 0 sind(a(1)) cosd(a(1))];
        RMy = [cosd(a(2)) 0 sind(a(2)); 0 1 0; -sind(a(2)) 0 cosd(a(2))];
        RMz = [cosd(a(3)) -sind(a(3)) 0; sind(a(3)) cosd(a(3)) 0; 0 0 1];
        RM  = RMx*RMy*RMz;
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
        % -- Also ensure eigenvectors point down.
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
%xmn = sum(xf(:,jc).*repmat((1-msft(jc)').^2,5,1),2)/sum((1-msft(jc)).^2);
xmn = sum(xf.*repmat((1-msft').^2,5,1),2)/sum((1-msft).^2);
X      = xmn([1 2 3; 2 4 5; 3 5 1]);
X(3,3) = -X(1,1)-X(2,2);
[V,S]  = eig(X);
S = flip(diag(S));
V = fliplr(V);
V(:,V(3,:)>0) = -V(:,V(3,:)>0);
PLNG = -asind(V(3,:))';
TRND = 90-atan2d(V(2,:),V(1,:))';
R = (S(1)-S(2))/(S(1)-S(3)); 

TRND(TRND<0) = TRND(TRND<0)+360;


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
%uR  = prctile(abs(Rs-R),prcntU(:));

bU = 50+[-pU'; pU']/2;
uR = reshape(prctile(Rs,bU(:)),[2 NU]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 2020-09-24 Addition: Contour Plots %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% -- The axes are stored in xyz1,xyz2,xyz3. Convert to plot coordinates
% -- THIS IS NOT CORRECT!!!!!
%xy1 = (xyz1(1:2,:).*repmat( sqrt(1-xyz1(3,:)),2,1))';
%xy2 = (xyz2(1:2,:).*repmat( sqrt(1-xyz2(3,:)),2,1))';
%xy3 = (xyz3(1:2,:).*repmat( sqrt(1-xyz3(3,:)),2,1))';

xyz1 = xyz1';
xyz2 = xyz2';
xyz3 = xyz3';
xyz1(xyz1(:,3)>0,:) = -xyz1(xyz1(:,3)>0,:);
xyz2(xyz2(:,3)>0,:) = -xyz2(xyz2(:,3)>0,:);
xyz3(xyz3(:,3)>0,:) = -xyz3(xyz3(:,3)>0,:);
azi1 = atan2d(xyz1(:,2),xyz1(:,1));
azi2 = atan2d(xyz2(:,2),xyz2(:,1));
azi3 = atan2d(xyz3(:,2),xyz3(:,1));
toa1 = acosd(xyz1(:,3));
toa2 = acosd(xyz2(:,3));
toa3 = acosd(xyz3(:,3));
r1   = sqrt(2)*cosd(toa1/2);
r2   = sqrt(2)*cosd(toa2/2);
r3   = sqrt(2)*cosd(toa3/2);
xy1  = repmat(r1,1,2).*[cosd(azi1) sind(azi1)]; 
xy2  = repmat(r2,1,2).*[cosd(azi2) sind(azi2)];
xy3  = repmat(r3,1,2).*[cosd(azi3) sind(azi3)];


N = 100;
xh = linspace(-1,1,N);
[xg,yg] = meshgrid(xh,xh);

X = xg*N;
Y = yg*N;
rad = sqrt(xg.^2+yg.^2);

%px1 = hist3(xy1,'Edges',{xh,xh})/nx;
%px2 = hist3(xy2,'Edges',{xh,xh})/nx;
%px3 = hist3(xy3,'Edges',{xh,xh})/nx;
%px1 = imgaussfilt(px1,1.0);
%px2 = imgaussfilt(px2,1.0);
%px3 = imgaussfilt(px3,1.0);


hst1 = hist3(xy1,'Edges',{xh,xh})/nx;
hst2 = hist3(xy2,'Edges',{xh,xh})/nx;
hst3 = hist3(xy3,'Edges',{xh,xh})/nx;

%H1 = imgaussfilt(hst1,1.1);
%H2 = imgaussfilt(hst2,1.5);
%H3 = imgaussfilt(hst3,1.5);

H1 = zeros(size(rad));
H2 = zeros(size(rad));
H3 = zeros(size(rad));
sgflt = N*0.025;
for ii = find(rad < 1)';
    
    jd = find( (sqrt(((X-X(ii)).^2+(Y-Y(ii)).^2)) < 3*sgflt).*(rad < 1));
    w = zeros(size(rad));
    w(jd) = exp(-0.5*sgflt^-2*( (X(jd)-X(ii)).^2+ (Y(jd)-Y(ii)).^2 ));
    w(jd) = w(jd)/sum(w(jd));

    H1(ii) = sum(w(jd).*hst1(jd));
    H2(ii) = sum(w(jd).*hst2(jd));
    H3(ii) = sum(w(jd).*hst3(jd));
    
end
%}
[px1,jsrt1] = sort(H1(:),'descend');
[px2,jsrt2] = sort(H2(:),'descend');
[px3,jsrt3] = sort(H3(:),'descend');
cx1 = cumsum(px1)/sum(px1);
cx2 = cumsum(px2)/sum(px3);
cx3 = cumsum(px3)/sum(px3);
cn1 = px1([find(cx1 > 0.10,1), find(cx1 > 0.50,1), find(cx1 > 0.90,1)]);
cn2 = px2([find(cx2 > 0.10,1), find(cx2 > 0.50,1), find(cx2 > 0.90,1)]);
cn3 = px3([find(cx3 > 0.10,1), find(cx3 > 0.50,1), find(cx3 > 0.90,1)]);
%cntr1 = nx*px1([find(cx1 > 0.10,1),find(cx1 > 0.34,1),find(cx1 > 0.6827,1)]);
%cntr2 = nx*px2([find(cx2 > 0.10,1),find(cx2 > 0.34,1),find(cx2 > 0.6827,1)]);
%cntr3 = nx*px3([find(cx3 > 0.10,1),find(cx3 > 0.34,1),find(cx3 > 0.6827,1)]);


figure(7002)
clf
hold on 
%imagesc(xh,xh,px1')
%colormap([ones(256,1),repmat((255:-1:0)',1,2)/255])

% -- Plot grid lines...
aa = [0:15:90];
ab = [0:30:150];
for ii = 1:length(aa)
    toa = aa(ii);
    azi = 1:360;
    r   = sqrt(2)*cosd((180-toa)/2);
    x   = r.*cosd(90-azi);
    y   = r.*sind(90-azi);
    plot(x,y,'Color',[0.7,0.7,0.7])
end
for ii = 1:length(ab)
    toa = -90:90;
    azi = ab(ii);
    r   = sqrt(2)*cosd((180-toa)/2);
    x   = r.*cosd(90-azi);
    y   = r.*sind(90-azi);
    plot(x,y,'Color',[0.7,0.7,0.7])
end
        
plot(cosd(0:360),sind(0:360),'k')
axis equal

[c1,h1] = contour(xg,yg,H1',cn1,'r');
[c2,h2] = contour(xg,yg,H2',cn2,'g');
[c3,h3] = contour(xg,yg,H3',cn3,'b');

% -- Add points for best solution
V(:,V(3,:)>0) = -V(:,V(3,:)>0);

azi = atan2d(V(2,:),V(1,:))';
toa = acosd(V(3,:))';
r   = sqrt(2)*cosd(toa/2);
xy  = repmat(r,1,2).*[cosd(azi) sind(azi)];
plot(xy(:,1),xy(:,2),'ko')


%imagesc(xh,xh,H1')

%[c1,h1] = contour(xg,yg,hst1,cntr1);
%cntrP = dscndP([find(cP > 0.10,1),find(cP > 0.34,1),find(cP > 0.6827,1)])








