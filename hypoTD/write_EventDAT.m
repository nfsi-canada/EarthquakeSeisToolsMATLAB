function write_EventDAT(file,evdir)
% function write_EventDAT(file,evdir)
% 
% 2020-10-13
% This function writes a hypoDD "event.dat"-style from an "evdir" cell
% array in Matlab.
%
%    INPUTS
% 
%       f == filename e.g. '...PATH/phase.dat'
%   evdir == event directory {T,LAT,LON,DEP,MAG,NAME}
%

Ne = size(evdir,1);
f = fopen(file,'w');
for ie = 1:Ne
    T   = evdir{ie,1};
    ymd = datestr(T,'yyyymmdd');
    hms = datestr(T,'hhMMssfff');
    hms = hms(1:end-1);
    lat = evdir{ie,2};
    lon = evdir{ie,3};
    dep = evdir{ie,4};
    mag = evdir{ie,5};
    fprintf(f,'%8s %8s %7.4f %8.4f %5.1f %3.1f 0.00 0.00 0.00 %d\n', ...
               ymd,hms,lat,lon,dep,mag,ie);
       
end
fclose(f);