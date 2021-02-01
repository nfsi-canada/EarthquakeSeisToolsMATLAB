% 2020-12-02
% This tests HypoAPP3D with a very simple velocity model...this
% had better work first before any more complicated models.

clear all

sigP = 0.0; % standard deviation (s) of pick noise
sigH = 5;    % standard deviation (km) of initial hypocenter errors
vpvs = 1.75;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Create velocity model
gx = [-100:10:100];
gy = [-80:10:80];
gz = [0:4:100];
[GX,GY,GZ] = ndgrid(gx,gy,gz);
V0 = 5 + sqrt(GZ)/8 + GX/500; 
V0 = V0 + min(1.2,1./( ((GX-0)/80).^2 + ((GY+100)/40).^2 + ((GZ-50)/40).^2));
V = griddedInterpolant({gx,gy,gz},V0);

xmin = min(gx);
ymin = min(gy);
xmax = max(gx);
ymax = max(gy);
xrng = xmax-xmin;
yrng = ymax-ymin;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Create stations
Ns = 10;
sta = zeros(Ns,3);
sta(:,1) = xrng*rand(Ns,1)+xmin;
sta(:,2) = yrng*rand(Ns,1)+ymin;
%sta(:,3) = rand(Ns,1).^3;
%sta(:,3) = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Create hypocenters
Ne = 10;
H0 = zeros(Ne,3);
H0(:,1) = (xrng*rand(Ne,1)+xmin)/2;
H0(:,2) = (yrng*rand(Ne,1)+ymin)/2;
H0(:,3) = 30*rand(Ne,1);

% -- Initial guesses should be "perturbed" a little
Hp = H0+sigH*randn(Ne,3);
Hp(Hp(:,3)<0,3) = 0.5;

T0 = repmat(datetime(2020,1,1,0,0,0),Ne,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Make picks
P0 = cell(Ne,1);
for ie = 1:Ne
    tpp = RayTrace3D(repmat(H0(ie,:),Ns,1),sta,V);
    P0{ie} = [sta   ones(Ns,1) ones(Ns,1) tpp
              sta 2*ones(Ns,1) ones(Ns,1) tpp*vpvs];    
end

% -- Add pick noise
Pp  = P0;
for ie = 1:Ne
    Pp{ie}(:,6) = Pp{ie}(:,6) + sigP*randn(2*Ns,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Estimate hypocenters
params.vpvs     = vpvs;
params.tol      = 1e-2;
params.picktol  = 6;
params.minZ     = 0.1;
params.maxZ     = 200;
params.Lcube    = 50;
params.NstepTol = 10;
params.minZ     = 0.1;
params.maxZ     = 250;
params.CI       = 0.9;
params.minpha   = 8;
params.max_iter = 200;
params.sclF     = 6;
%params.steplen = 2;

S1 = cell(Ne,1);
for ie = 1:Ne
    %[H1(ie,:),T1(ie),P1{ie},S1{ie}] = hypoAPP3Dgrid(Hp(ie,:),T0(ie),Pp{ie},V,params);
    [H1(ie,:),T1(ie),P1{ie},S1{ie}] = hypoAPP3D(Hp(ie,:),T0(ie),Pp{ie},V,params);
end

[sqrt(sum((H0-Hp).^2,2)) sqrt(sum((H0-H1).^2,2))]
