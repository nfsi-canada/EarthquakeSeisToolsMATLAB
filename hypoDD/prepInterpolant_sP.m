function [FZsp,FTsp] = prepInterpolant_sP(model,vpvs,Xmax,Zmax)
% function [FZsp,FTsp] = prepInterpolant_sP(model,vpvs,Xmax,Zmax)
%
% 2021-01-16
% This function creates an interpolant that returns hypocenter depth
% given a epicentral distance and sP-minus-P time. Note that this does
% not take topography / stations elevation into account, it assumes the 
% bounce-point and station are at "zero".
%
% INPUTS
%
%   model = [Z(km), V(km s-1)] N_layer x 2 matrix
%    vpvs = (uniform) Vp/Vs ratio
%    Xmax = maximum epicentral distance
%    Zmax = maximum event depth
%
% OUTPUTS
%
%     FZsp = interpolant that returns depth given epicentral distance
%            and sP time. [Z] = FZsp(X,T (sP-P))
%     FTsp = interpolant that return sP-P time given distance and depth
%            [TsP] = FTsp(X,Z)

Nx = 300;
Nz = 300;
xrng = [0; exp(linspace(0.1,log(Xmax),Nx-1)')-1];
zrng = exp(linspace(0.1,log(Zmax),Nz)')-1;    



% -- Prep for sP depth constraints Bounce P-times for full distance range
tb = RayTrace(1e-6,xrng,model);

xS  = zeros(Nx,Nx);
TsP = zeros(Nx,Nz);


for iz = 1:Nz
    
    % -- Get direct P travel times
    tP0 = RayTrace(zrng(iz),xrng,model);
   
    % -- For each X, need to try all the potential bounce points
    % -- and find point with minimum sP time. For each TOTAL ditance, 
    % -- interpolate Nx bounce-point distances. 
    % -- xS == source to bounce-point distance 
    % -- xP == bounce-point to receiver distance
    for ix = 1:Nx
        xS(:,ix) = linspace(0,xrng(ix),Nx);       
    end
    xP = repmat(xrng',Nx,1)-xS;

    tS_up  = vpvs*interp1(xrng,tP0,xS);
    tP_bnc = interp1(xrng,tb,xP);
    
    % -- Ideally I would find subsample minimum here
    TsP(:,iz) = min(tS_up + tP_bnc)'-tP0;
  
end

% -- For some very short distances sP-P is negligible and computing
% -- errors will have made it negative, correct this
TsP(TsP<0) = 0;

% -- Create interpolant
X = repmat(xrng,Nz,1);
Z = repmat(zrng',Nx,1);
FZsp = scatteredInterpolant(X(:),TsP(:),Z(:));
FTsp = scatteredInterpolant(X(:),Z(:),TsP(:));

