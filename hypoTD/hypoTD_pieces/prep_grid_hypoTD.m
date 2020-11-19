% 2020-05-06
% This piece of hypoTD sets up a grid of travel-times, take-off angles,
% and local velocities to be interpolated


% -- Make tables (distance x depth) of travel times, ray angles 
d0 = distance(repmat(hyp0(:,1),1,Ns),repmat(hyp0(:,2),1,Ns), ...
     repmat(stadir(:,1)',Ne,1),repmat(stadir(:,2)',Ne,1),[6371,0]);
maxD = max(d0(:)); 

drng = linspace(0.001,maxD,500)';
zrng = linspace(minZ,maxZ,200)';
nx   = length(drng);
nz   = length(zrng);
[Zrng,Drng] = meshgrid(zrng,drng);

tp  = zeros(nx,nz);
tht = zeros(nx,nz);
vz  = zeros(nz,1);
for iz = 1:nz
    [tp(:,iz),tht(:,iz),~,~,vz(iz)] = RayTrace(zrng(iz),drng,model);
end


% -- DO I NEED TO USE A MAP PROJECTION?? I THINK SO
lat0 = min(hyp0(:,1))-0.5;
lat1 = max(hyp0(:,1))+0.5;
lon0 = min(hyp0(:,2))-0.5;
lon1 = max(hyp0(:,2))+0.5;

map = m_proj('UTM', 'lat',[lat0,lat1],'long',[lon0,lon1],'ell','wgs84','rect','on');


% -- Prep for sP depth constraints Bounce P-times for full distance range
tBounce = RayTrace(1e-6,drng,model);

% -- For each event, store indices of all its phases,
% -- as well as their demeaned travel times
tph0 = zeros(Np,1);
jphE = cell(Ne,1);
for ie = 1:Ne
    jphE{ie} = find(phadir(:,1)==ie);
    tph0(jphE{ie}) = phadir(jphE{ie},4)-mean(phadir(jphE{ie},4));
end
    
