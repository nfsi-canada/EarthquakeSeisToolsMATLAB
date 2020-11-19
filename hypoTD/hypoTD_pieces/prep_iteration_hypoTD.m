% 2020-05-05
% This piece of hypoTD prepares to take invert triple-differences
% by computing inter-event distances, ray-angles to each station, etc.
%
%


% -- Convert lat/lon to map coorditnates (in km)
% -- Need to check for events outside map bounds...

[ex,ey] = m_ll2xy(hyp(:,2),hyp(:,1));

jnan = find(isnan(ex));
hypNaN = hyp(jnan,1:2);
hypNaN(hypNaN(:,1) < lat0,1) = lat0;
hypNaN(hypNaN(:,1) > lat1,1) = lat1;
hypNaN(hypNaN(:,2) < lon0,2) = lon0;
hypNaN(hypNaN(:,2) > lon1,2) = lon1;
[ex(jnan),ey(jnan)] = m_ll2xy(hypNaN(:,2),hypNaN(:,1));

mx = mean(ex);
my = mean(ey);
ex = (ex-mx)/1e3;
ey = (ey-my)/1e3;
ez = hyp(:,3);

%if iter==1
%    disp('HypoTD Starting from Mean!!!')
%    ex(:) = mean(ex);
%    ey(:) = mean(ey);
%    ez(:) = mean(ez);
%end
x0 = [ex; ey; ez];

% -- Epicentral distance, azimuth, from (Event x Station)
[dES,aES] = distance(repmat(hyp(:,1),1,Ns),repmat(hyp(:,2),1,Ns), ...
                     repmat(stadir(:,1)',Ne,1),repmat(stadir(:,2)',Ne,1),[6371,0]);
zES = repmat(hyp(:,3),1,Ns);     

% -- Predict times, angles, by interpolation of previous table
tppES = interp2(Zrng,Drng,tp, zES,dES); 
thtES = interp2(Zrng,Drng,tht,zES,dES);

% -- Sines,cosines (note that DOWN is the positive vertical direction)
cH = cosd(90-aES);
sH = sind(90-aES);
cV = cosd(180-thtES);
sV = sind(180-thtES);

% -- Find comparisons of events within max. distance
% -- put in separate matrices specific to this iteration

dCC = sqrt(distance(hyp(mttcc(:,1),1),hyp(mttcc(:,1),2), ...
                    hyp(mttcc(:,2),1),hyp(mttcc(:,2),2),[6371 0]).^2 ...
           +(hyp(mttcc(:,1),3)-hyp(mttcc(:,2),3)).^2);
dCT = sqrt(distance(hyp(mttct(:,1),1),hyp(mttct(:,1),2), ...
               hyp(mttct(:,2),1),hyp(mttct(:,2),2),[6371 0]) ...
           +(hyp(mttct(:,1),3)-hyp(mttct(:,2),3)).^2);         
jcc = find(dCC < mxDCC(iter));
jct = find(dCT < mxDCT(iter));
MCC = mttcc(jcc,:);
MCT = mttct(jct,:);
dCC = dCC(jcc);
dCT = dCT(jct);

Ncc = size(MCC,1);
Nct = size(MCT,1);
NsP = size(sP,1);

