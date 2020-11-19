% 2019-10-31
% This piece randomly generates the locations and MTs for a synthetic relMT
% data set. Important outputs are ex,ey,ez, M, and 'm' which contains the 6
% independent elements of each event in one big matrix. 

%if ~exist('sdrs')


latCEN = 44.639;
lonCEN = -63.570;
TCEN   = datetime(2020,01,01,00,00,00);

% -- Event locations
ex = sqrt(EvSp)*randn(ne,1);
ey = sqrt(EvSp)*randn(ne,1);
ez = sqrt(EvSp)*randn(ne,1) + EvDp;

% -- Could add 'Noise' factor to these locations hypocenter estimates
exE = ex + sqrt(EvHE)*randn(ne,1) + sqrt(EvHS)*randn;
eyE = ey + sqrt(EvHE)*randn(ne,1) + sqrt(EvHS)*randn;
ezE = ez + sqrt(EvHE)*randn(ne,1) + sqrt(EvHS)*randn;

[latE,lonE] = reckon(latCEN,lonCEN,sqrt(exE.^2+eyE.^2),90-atan2d(eyE,exE),[6371 0]);
[latTRUE,lonTRUE] = reckon(latCEN,lonCEN,sqrt(ex.^2+ey.^2),90-atan2d(ey,ex),[6371 0]);

hyp0    = [latE,lonE,ezE];
hypTRUE = [latTRUE,lonTRUE,ez];

T0 = repmat(TCEN,ne,1)+days(0:ne-1)';
%{

    
% -- Event magnitudes
%Mw = 3+0.2*randn(ne,1);
%Mo = mom2mag(Mw,1);
%MoSclFac = 1000/mean(Mo);
%MoNrm = MoSclFac*Mo;

% -- EQ fault parameters   
str0 = 180;
dip0 = 45;
rke0 = 0;
slp0 = 0;
sdrs = [str0*ones(ne,1), dip0*ones(ne,1), rke0*ones(ne,1), slp0*ones(ne,1)]; 
sdrs(:,1) = sdrs(:,1) + rngS * (2*(rand(ne,1)-0.5));
sdrs(:,2) = sdrs(:,2) + rngD * (2*(rand(ne,1)-0.5));
sdrs(:,3) = sdrs(:,3) + rngR * (2*(rand(ne,1)-0.5));
sdrs(:,4) = sdrs(:,4) + rngA * (2*(rand(ne,1)-0.5));



% -- Convert fault parameters + magnitudes to MTs
M0 = cell(ne,1);
m  = zeros(6,ne);
for ie = 1:ne
    M0{ie}   = MoNrm(ie)*sdrs2mij(sdrs(ie,:));
    m(:,ie) = M0{ie}([1 5 9 2 3 6]');
end

% -- All events are used in the synthetic
nue = ne;
ue  = [1:ne]';
uMo = MoNrm;

nref = length(jref);
jrf0 = jref;
gmref = m(:,jref);
nrfT  = length(jref);
x06 = m(:);
x05 = m([1:2,4:end],:);
x05 = x05(:);
%}