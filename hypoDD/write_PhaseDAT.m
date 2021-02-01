function write_PhaseDAT(file,evdir,stadir,phadir,je)
% function write_PhaseDAT(file,evdir,stadir,phadir,je)
% 
% 2020-10-13
% This function writes a hypoDD "phase.dat"-style from event, station,
% and phase directories.
%
%    INPUTS
% 
%       f == filename e.g. '...PATH/phase.dat'
%   evdir == event directory {T,LAT,LON,DEP,MAG,NAME}
%  stadir == station directory {NAME,LAT,LON,ELV}
%  phadir == [EV#,STA#,PHA#,T], PHA# = 1/2 for P/S
%      je == (optional) list of event indices to include. Must match 
%            size of evdir. phadir(:,1) must match these indices.
%

Ne = size(evdir,1);
if nargin < 5;
    je = 1:Ne;
end

f = fopen(file,'w');
PS = {'P','S'};
for ie = 1:Ne
    T   = evdir{ie,1};
    lat = evdir{ie,2};
    lon = evdir{ie,3};
    dep = evdir{ie,4};
    mag = evdir{ie,5};
    fprintf(f,'# %4d %02d %02d %02d %02d %02.0f %7.4f %8.4f %5.1f %3.1f 0 0 0 %d\n', ...
            T.Year,T.Month,T.Day,T.Hour,T.Minute,T.Second,lat,lon,dep,mag,je(ie));
    
    phaE = phadir(phadir(:,1)==je(ie),:);
    for ii = 1:size(phaE,1)
        sta = stadir{phaE(ii,2),1};
        T   = phaE(ii,4); 
        W   = phaE(ii,5);
        ph  = PS{phaE(ii,3)};
        fprintf(f,'%7s %6.2f %5.3f %1s\n',sta,T,W,ph);
    end
end
fclose(f);