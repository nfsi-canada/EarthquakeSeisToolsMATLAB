% 2019-10-31
% This piece generates a complete set of possible P- and S-  phases given
% sets of events/stations from createEvents/createStations. Important outputs
% include 'phadir', which contains the arrival times, as well as:
%  azi = 'azimuth' of ray in degrees
%  toa = 'take-off angle' of ray in degrees, where 0=up, 90=horiz., 180=down
%  rdp = 'ray-path distance' in km
%   ra = 'receiver angle', andlge of incidence at receiver, again 0=up
%
%


% -- I should do this all twice...once for actual locations to make waveforms,
% -- and again with 'estimated' locations, which have noise added


% -- Now need azimuth, take-off angle, ray-path distance, for each event
% -- to each station. Plus I need the receiver angle
azi = zeros(ne,ns);
toa = zeros(ne,ns);
rpd = zeros(ne,ns);
tP  = zeros(ne,ns);
ra  = zeros(ne,ns);
for ie = 1:ne
    aziE = 90-atan2d(sy-ey(ie),sx-ex(ie));
    aziE(aziE < 0)   = aziE(aziE < 0) + 360;
    aziE(aziE > 360) = aziE(aziE > 360) - 360;
    azi(ie,:) = aziE;
    [tP(ie,:),toa(ie,:),rpd(ie,:),ra(ie,:)] = RayTrace(ez(ie),sqrt((ex(ie)-sx).^2+(ey(ie)-sy).^2),model);
end
tS = sqrt(3)*tP;
col = @(X) X(:);

% -- Construct phadir [EV,STA,PHA,T], phaprms [azi,toa,rpd]
jE = col(repmat(1:ne,2*ns,1));
jS = repmat(col(repmat(1:ns,2,1)),ne,1);
jP = repmat([1;2],ne*ns,1);
tt = col([col(tP')'; col(tS')']);

phadir  = [jE,jS,jP,tt];

% -- Randomly cull some?
[phadir,jkeep] = cullPhadirHypoTD(phadir,NF0);
np  = size(phadir,1);


%{
azi = col(repmat(col(azi')',2,1));
toa = col(repmat(col(toa')',2,1));
rpd = col(repmat(col(rpd')',2,1));
ra  = col(repmat(col(ra')',2,1));

% -- True ray angles/distance in 'phaprms0' 
phaprms0 = [azi,toa,rpd,ra];
phaprms0 = phaprms0(jkeep,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Repeat for ESTIMATED hypocenters

% -- Store Estimated ray angles/distance in 'phaprms' rather than phaprms0
% -- Don't need to re-do cullPhadir. 

azi = zeros(ne,ns);
toa = zeros(ne,ns);
rpd = zeros(ne,ns);
tP  = zeros(ne,ns);
ra  = zeros(ne,ns);
for ie = 1:ne
    aziE = 90-atan2d(sy-eyE(ie),sx-exE(ie));
    aziE(aziE < 0)   = aziE(aziE < 0) + 360;
    aziE(aziE > 360) = aziE(aziE > 360) - 360;
    azi(ie,:) = aziE;
    [tP(ie,:),toa(ie,:),rpd(ie,:),ra(ie,:)] = RayTrace(ezE(ie),sqrt((exE(ie)-sx).^2+(eyE(ie)-sy).^2),'GSC');
end

azi = col(repmat(col(azi')',2,1));
toa = col(repmat(col(toa')',2,1));
rpd = col(repmat(col(rpd')',2,1));
ra  = col(repmat(col(ra')',2,1));

phaprms = [azi,toa,rpd,ra];
phaprms = phaprms(jkeep,:);
%}