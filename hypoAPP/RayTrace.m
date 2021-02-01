function [tpp,theta,dray,thetaR,v_hyp] = RayTrace(hdepth,xep,model)
% function [tpp,theta,dray,thetaR,v_hyp] = RayTrace(hdepth,xep,model)
%
% 2020-03-13
% This function finds the travel time for a given depth, distance, 
% in a gradient velocity model. It also returns the takeoff angle in
% degrees, measured from vertical (up). It can take in multiple xep 
% (station distances) but only for one EQ (one hdepth). 
%
%     INPUTS
%
%      hdepth == hypocentral depth in km
%         xep == epicentral distance in km
%       model == name of velocity model e.g. 'GSC','ALASKA','ma2011'
%                OR a table of [Depth (km), Velocity (km s-1)]
%                this takes in a gradient model, not a layer cake
%
%     OUTPUTS
%
%         tpp == predicted P travel time in s
%       theta == the take-off angle, 0 being vertically upward
%        dray == ray-path distance in km
%      thetaR == angle of incidence at receiver, 0 being vertically upward
%       v_hyp == P-wave velocity at the hypocenter 
%                (TRUE VELOCITY, no earth-flattening transformation)

if ischar(model)
    switch model
    case 'ma2011'
    DV = [  0.0  6.25
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
    case 'GSC'
        DV = [0.0  5.00
          1.0  5.50
          6.0  6.35
         30.0  6.90
         45.0  7.42
         50.0  7.77
         60.0  7.78
        100.0  8.00
        600.0 12.00];
    case 'ALASKA'
        DV = [0.00  5.8
          24.3  6.0
          24.4  7.3
          40.1  7.5
          40.2  7.8
          75.9  8.0
          76.0  8.27
         300.0  8.31
         301.0  10.3
         524.0  10.5       
         525.00 12.6
        1000.00 13.0];
    case 'customLSL'
        DV = [0  6.07
              8  6.45
             17  6.60
             24  6.69
             30  6.83
             36  7.04
             38  7.25
             40  7.80
             45  7.87
             60  7.99
            200  8.25
            600  9.50]; 
    end
else
    DV = model;
end

depth = DV(:,1);
vp    = DV(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a   = 6371; % km
rf  = a-depth;
zf  = -a*log(rf/a);
vf  = (a./rf).*vp;
hzf = -a*log((a-hdepth)/a);


id = find(zf>hzf,1);


% -- Find last velocity w/ linear gradient, relative EQ depth in layer
hv = vf(id-1)+(vf(id)-vf(id-1))*(hzf-zf(id-1))/(zf(id)-zf(id-1));
zf = [zf(1:id-1); hzf; zf(id:end)];
vf = [vf(1:id-1);  hv; vf(id:end)];    
 

lnumUP = id-1;
lnumT  = length(zf)-1;

% -- 
np   = 200;
tht = linspace(0.1,89.9,np)'; 
p   = sind(tht)/hv;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Up-going rays
X1 = zeros(np,1);
T1 = zeros(np,1);
D1 = zeros(np,1);

for ip = 1:np
    x = 0;
    t = 0; 
	d = 0;
    for iL = 1:lnumUP
        h  = zf(iL+1)-zf(iL);
        vt = vf(iL);
        vb = vf(iL+1);
        [dx,dt,irtr] = layerxt(p(ip),h,vt,vb);
    
        x = x + dx;
        t = t + dt; 
		d = d + sqrt(dx^2+h^2);
    end
    X1(ip) = x;
    T1(ip) = t;
	D1(ip) = d;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Interpolations

% -- Method for down-going rays isn't very stable, so just extrapolate 
% -- from upgoing rays
%tpp    = interp1(X,T,xep,'pchip','extrap');
%dray   = interp1(X,D,xep,'pchip','extrap');
%theta  = interp1(X,tht,xep,'linear','extrap');

%theta(theta>179) = 179 + (theta(theta>179)-179)/(max((theta(theta>179)-179))+0.1);
%thetaR = asind((vf(1)/hv)*sind(theta));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Down-going rays

X2 = zeros(np,1);
T2 = zeros(np,1);
D2 = zeros(np,1);
for ip = 1:np 
    x = 0;
    t = 0; 
	d = 0;
    for iL = id:lnumT 
        h  = zf(iL+1)-zf(iL);      
        vt = vf(iL);
        vb = vf(iL+1);
        [dx,dt,irtr] = layerxt(p(np+1-ip),h,vt,vb);

        x = x + dx;
        t = t + dt; 
		d = d + sqrt(dx^2+h^2);			

        if irtr == 2 
            break
        end       
    end
    %if irtr~=2
    %    keyboard
    %    break
    %end
    X2(ip) = x;
    T2(ip) = t;
	D2(ip) = d;
end

X2 = 2*X2+flipud(X1);
T2 = 2*T2+flipud(T1);
D2 = 2*D2+flipud(D1);

X      = [X1; X2];
T      = [T1; T2];
D      = [D1; D2];
tht    = [tht; 90+tht];


% -- Check for decreases in X, create a smoothly increasing X vs. tht plot
% -- This method is imperfect because we don't necessarilly know
% -- which ray of a triplication is faster (deep vs. shallow), but this
% -- code assumes the shallower one is faster. This routine will return
% -- a shallower take-off angle then the "bending" routine in RayTrace1D
% -- in some (relattively rare) cases
while 1 
    jn = find(diff(X)<0)+1;
    if ~length(jn)
        break
    end
    X(jn) = [];
    T(jn) = [];
    D(jn) = [];
    tht(jn) = [];
end


% -- Should find minimum arrival time for a given distance
tpp   = interp1(X,T,xep,'pchip','extrap');
dray  = interp1(X,D,xep,'pchip','extrap');
theta = interp1(X,tht,xep,'linear','extrap');
theta(theta>179) = 179 + (theta(theta>179)-179)/(max((theta(theta>179)-179))+0.1);
thetaR = asind( (vf(1)/hv)*sind(theta) );

% -- P-wave veloctiy in km s-1, no earth-flattening transformation
v_hyp = interp1(depth,vp,hdepth);

% -- DO I NEED TO ACCOUNT FOR EARTH-FLATTENING TRANSFORMATION 
% -- IN RAY-ANGLE PREDICTION??
