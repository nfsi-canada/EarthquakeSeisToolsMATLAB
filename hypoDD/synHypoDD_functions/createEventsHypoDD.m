function [hyp0,hyp1] = createEventsHypoDD(Ne,latC,lonC,depC,scl,errR,errS)
% function [hyp0,hyp1] = createEventsHypoDD(Ne,latC,lonC,depC,scl,errR,errS)
%
% 2021-01-11
% This piece randomly generates the locations and MTs for a synthetic relMT
% data set. Important outputs are ex,ey,ez, M, and 'm' which contains the 6
% independent elements of each event in one big matrix.
%
% INPUTS
% 
%     Ne  == number of events
%    latC == centroid latitude
%    lonC == centroid longitude
%    depC == centroid depth (km)
%     scl == scale of cluster (km), standard deviation of event coordinates
%    errR == standard deviation of errors to individual events (km)
%    errS == absolute error in cluster centroid shift (km)
%
% OUTPUTS
%
%    hyp0 == true event locations [lat,lon,dep]
%    hyp1 == initial/perturbed event locations [lat,lon,dep]


% -- Event locations
ex = sqrt(scl)*randn(Ne,1);
ey = sqrt(scl)*randn(Ne,1);
ez = sqrt(scl)*randn(Ne,1) + depC;
ez(ez < 0) = -ez(ez < 0);

% -- Could add 'Noise' factor to these locations hypocenter estimates
shftC = randn(3,1);
shftC = errS*shftC/norm(shftC);
ex1 = ex + sqrt(errR)*randn(Ne,1) + shftC(1);
ey1 = ey + sqrt(errR)*randn(Ne,1) + shftC(2);
ez1 = ez + sqrt(errR)*randn(Ne,1) + shftC(3);
ez1(ez1 < 0) = -ez1(ez1 < 0); 

[lat1,lon1] = reckon(latC,lonC,sqrt(ex1.^2+ey1.^2),90-atan2d(ey1,ex1),[6371 0]);
[lat0,lon0] = reckon(latC,lonC,sqrt(ex.^2+ey.^2),90-atan2d(ey,ex),[6371 0]);

hyp0 = [lat0,lon0,ez];
hyp1 = [lat1,lon1,ez1];


