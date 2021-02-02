function [AR,bR,NR] = build_AFM(PVH,Qp,Qs,vpvs)
% function [AR,bR,NR] = build_AFM(PVH,Qp,Qs,vpvs)
%
% 2020-11-19
% This function makes a linear system out of P-wave first-motion data,
% where the solution vector would contain unit moment tensor elements.
%
% For vertical-only stations, it "corrects" amplitudes on the vertical  
% component to P,SV amplitudes, which can be used if the 
% Kisslinger [1981,BSSA] (P/SV)_z method is used
%
%  INPUTS
%
%    PVH = [AZI,TOA,TP,TS,FRQ,Ap,Av,Ah]
%            AZI == azimuth, degrees east of north
%            TOA == take-off angle in degrees
%                     0 == up 
%                    90 == horizontal
%                   180 == down
%            TP == P travel time (s)
%            TS == S travel time (s)
%           FRQ == Nominal frequency of bandpass amplitudes were 
%                  measured at (I recommend using the highpass frq.)
%            Ap == P absolute amplitude  
%            Av == SV absolute amplitude
%            Ah == SH absolute amplitude (or zero for single-component)
%     Qp = P-wave attenuation (unitless) typically order 1000 
%     Qs = S-wave attenuation (unitless) typically ~0.5*Qp
%   vpvs = average P-wave/S-wave velocity ratio
%   
%                  
%            
%
%  OUTPUTS
%
%   AR = an NR x 6 matrix of coefficients corresponding to the 
%        moment tensor elements [M11 M22 M33 M12 M13 M23]
%   bR = NR x 1 right-hand side vector of +1 or -1
%   NR = number of amplitude ratios able to be formed from PVH

% -- Number of independent MT elements [M11 M22 M33 M12 M13 M23]
nx = 6;

% -- Separate columns of PVH for ease of coding
azi = PVH(:,1);
toa = PVH(:,2);
TP  = PVH(:,3);
TS  = PVH(:,4);
FRQ = PVH(:,5);
Ap  = PVH(:,6);
Av  = PVH(:,7);
Ah  = PVH(:,8);

% -- Identify single-component (vertical) stations
jjV = find(Ah==0);
jj3 = find(Ah);
Ntmp = size(PVH,1);
NR = 3*length(jj3) + length(jjV);

% -- Direction cosines (E,N,U)...
gx = cosd(90-azi).*sind(toa);
gy = sind(90-azi).*sind(toa);
gz = cosd(toa);

% -- Vertical, SV 
zr = zeros(Ntmp,1);

gP  = [gx.^2, gy.^2, gz.^2, 2*gx.*gy, 2*gx.*gz, 2*gy.*gz];
gSx = [gx, zr, zr, gy, gz, zr] - repmat(gx,1,6).*gP;
gSy = [zr, gy, zr, gx, zr, gz] - repmat(gy,1,6).*gP;
gSz = [zr, zr, gz, zr, gx, gy] - repmat(gz,1,6).*gP;

gSH = zeros(Ntmp,nx);
gSV = zeros(Ntmp,nx);
for ii = 1:nx
    gR        = cosd(90-azi).*gSx(:,ii) + sind(90-azi).*gSy(:,ii);
    gSH(:,ii) = sind(90-azi).*gSx(:,ii) - cosd(90-azi).*gSy(:,ii);
    gSV(:,ii) = sind(toa).*gSz(:,ii) - cosd(toa).*gR;
end


% -- Weight by inverse square root of P-wave travel time. 
% -- Down-weight single component measurements by a factor of 3
% -- Upweight P-SH (the most stable ratio) by a factor of 3
wghtR      = TP.^(-1/2); 
wghtR      = [wghtR(jj3); 3*wghtR(jj3); wghtR(jj3); wghtR(jjV)/3];
wghtR      = wghtR/mean(wghtR); 

% -- 2020-07-13  MISSING (VP/VS)^2 ratio term????
% -- Or should it be cubed??
% -- Should I account for differential attenuation???
% -- Take highpass frequency used in record_P_SV_SH_amplitudes
NcyP  = FRQ.*TP;
NcyS  = FRQ.*TS;
dEP   = (1-2*pi/Qp).^NcyP/2; % Because amplitude prop. to energy^0.5
dES   = (1-2*pi/Qs).^NcyS/2;
PSfac = vpvs^3 * (dES./dEP);

% -- Build matrix 
AR = [gP(jj3,:),  gSV(jj3,:) 
      gP(jj3,:),  gSH(jj3,:)
      gSV(jj3,:), gSH(jj3,:)
      gP(jjV,:),  gSV(jjV,:)];  
  
bR = [(Ap(jj3)./Av(jj3)).*PSfac(jj3)
      (Ap(jj3)./Ah(jj3)).*PSfac(jj3)
      (Av(jj3)./Ah(jj3))
      (Ap(jjV)./Av(jjV)).*PSfac(jjV)];

% -- Make sure the predicted ratio is < 1  
jb1      = find(bR>1);
jb0      = find(bR<=1);
bR(jb1)  = bR(jb1).^-1;  
AR(jb1,:) = AR(jb1,[nx+1:2*nx,1:nx]);   

% -- Apply weights
nrmAb = vecnorm([AR bR]')';
AR = AR.*repmat(wghtR./nrmAb,1,2*nx);
bR = bR.*(wghtR./nrmAb);

