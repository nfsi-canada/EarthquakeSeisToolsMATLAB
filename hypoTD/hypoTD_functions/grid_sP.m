function TsP = sP_minus_P(x,z,model)
% function TsP = sP_minus_P(x,z,model)
%
% This function estimates the relative arrival time of the sP phase
% (sP - P) in seconds.
%
% I want to quickly compute sP times for range of depths, distances...
% I already have a grid of tp(x,z)

nz = length(z);
nx = length(x);

TsP = zeros(nx,nz);
tp2 = RayTrace(1e-6,flipud(drng),model);

for iz = 1:nz
   
   
    [tp1,theta1,dray,thetaRS] = RayTrace(z(iz),xrng,model);
    ts1 = tp1*sqrt(3);
    tp  = tp1(end);
    
    % -- Constant ray parameter at surface p = uP*sind(thetaP) = uS*sind(thetaS) 
    % -- For an S-to-P conversion... sind(thetaP) = sqrt(3)*sind(thetaS)
    thetaRP = asind( min(sind(thetaRS)*sqrt(3),1) );
    
    
    [tp2,theta2,dray,thetaR] = RayTrace(1e-6,flipud(xrng),model);

    tsP = ts1+tp2 - tp;
    
    [TsP(iz),jmin] = min(tsP);
    
end