function CC = DCpairCC(sdr1,sdr2,Np)
% function CC = DCpairCC(sdr1,sdr2,Np)
%
% 2020-07-06
% This function computes a "cross-correlation" coefficient between pairs
% of double-couple focal mechanisms. The CC is computed from their P-wave
% radiation patterns, after evenly sampling the focal sphere. 
%
%   INPUTS
%
%      sdr1 = [strike,dip,rake] 
%      sdr2 = [strike,dip,rake] (must be same height as sdr1)
%
%  OUTPUTS
%
%       CC = P-wave radiation cross-correlation coefficient

N = size(sdr1,1);

if nargin < 3
    Np = 2048;
end

j6 = [1 5 9 2 3 6]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Evenly sample focal sphere
% -- using fibonacci method

if mod(Np,2)
    Np = Np+1;
end

GoldenRatio = (1 + sqrt(5))/2;

ii  = [1:2*Np]'-0.5;
tht = acos(1-2*ii/(2*Np));
phi = 2*pi*ii/GoldenRatio;
phi = atan2(sin(phi),cos(phi));
jj  = find(tht < pi/2);

tht(jj) = [];
phi(jj) = [];

Np = length(jj);
gx = sin(tht).*cos(phi);
gy = sin(tht).*sin(phi);
gz = cos(tht);

g = [gx.^2, gy.^2, gz.^2, 2*gx.*gy, 2*gx.*gz, 2*gy.*gz];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = zeros(N,1);
for ii = 1:N
    
    M1 = sdrs2mij([sdr1(ii,:) 0]);
    M2 = sdrs2mij([sdr2(ii,:) 0]);
    x1 = M1(j6);
    x2 = M2(j6);
    
    CC(ii) = CCcoef(g*x1,g*x2);
end
  



