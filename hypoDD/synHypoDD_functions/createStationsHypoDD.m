function stadir = createStationsHypoDD(Ns,No,model,latC,lonC,depC)
% function stadir = createStationsHypoDD(Ns,No,model,latC,lonC,depC)
%
% 2021-01-11
% This piece randomly generates station in the 'octants' indicated by 'no'
% and computes their XY / EN coordinates, given a 3D velocity model
%
% INPUTS
%
%     Ns = number of stations
%     No = number octants stations exist in, i.e. No==4 for anywhere in 
%          upper hemisphere
%  model = strucutre with fields "mapcmd" and "V"   
%   latC = earthquake centroid latitude
%   lonC = earthquake centroid longitude
%   depC = earthquake centroid depth (km)
% 
% OUTPUTS 
%
%   stadir = Ns x 3 array [lat,lon,depS(km)]


% -- Station angles (azimuth, take-off angle) from event centroid
if No < 5
    aziS = 90*min(No,4)*rand(Ns,1);
    toaS = 90*sqrt(rand(Ns,1));
else
   aziS = zeros(Ns,1);
   toaS = zeros(Ns,1);
   for iss = 1:ns
      while 1
          azi = 360*rand;
          toa = 180*rand;
          jo  = 4*(toa > 90)+fix(azi/90) + 1;
          if jo <= No
              aziS(iss) = azi;
              toaS(iss) = toa;   
              break
          end
      end
   end
end

[cx,cy] = m_ll2xy(lonC,latC);
cx = cx/1000;
cy = cy/1000;
cz = depC;

% -- Do Ray-Tracing
dG       = [0, 10.^[0:0.1:4]]';
[~,toaG] = RayTrace(cz,dG,model);
dS       = interp1(toaG,dG,toaS,'pchip');

sx = dS.*cosd(90-aziS)+cx;
sy = dS.*sind(90-aziS)+cy;
[lonS,latS] = m_xy2ll(sx*1000,sy*1000);
stadir = [latS,lonS,-(rand(Ns,1)-0.1)];