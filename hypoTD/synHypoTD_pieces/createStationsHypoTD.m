% 2019-10-31
% This piece randomly generates station in the 'octants' indicated by 'no'
% and computes their XY / EN coordinates. 

% -- Station angles (azimuth, take-off angle) from event centroid
if no < 5
    aziS = 90*min(no,4)*rand(ns,1);
    toaS = 90*sqrt(rand(ns,1));
    %toaS = 89*ones(ns,1)
else
   aziS = zeros(ns,1);
   toaS = zeros(ns,1);
   for iss = 1:ns
      while 1
          azi = 360*rand;
          toa = 180*rand;
          jo  = 4*(toa > 90)+fix(azi/90) + 1;
          if jo <= no
              aziS(iss) = azi;
              toaS(iss) = toa;   
              break
          end
      end
   end
end

% -- Station distance from event centroid
% -- Need ray tracing code for super-simple velocity model
DepiGrd = [0, 10.^[0:0.1:4]]';
[~,toaGrd,rpdGrd] = RayTrace(EvDp,DepiGrd,model);
DepiS = interp1(toaGrd,DepiGrd,toaS,'pchip');

% -- Station coordinates
sx = DepiS.*cosd(90-aziS);
sy = DepiS.*sind(90-aziS);


[latS,lonS] = reckon(latCEN,lonCEN,DepiS,aziS,[6371 0]);
staLL = [latS,lonS];