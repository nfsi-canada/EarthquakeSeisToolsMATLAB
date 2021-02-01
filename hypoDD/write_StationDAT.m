function write_StationDAT(file,stadir)
% function write_StationDAT(file,stadir)
% 
% 2020-12-04
% This function writes a hypoDD "station.dat"-style from a cell array
% station directory
%
%    INPUTS
% 
%       f == filename e.g. '...PATH/phase.dat'
%  stadir == station directory {NAME,LAT,LON,ELV}
%

[NS,Nc] = size(stadir);

f = fopen(file,'w');

if Nc == 3
    for iss = 1:NS
        sta = pad(stadir{iss,1},7);
        lat = pad(num2str(stadir{iss,2},'%7.4f'),8,'left');
        lon = pad(num2str(round(stadir{iss,3}),'%8.4f'),9,'left');
        fprintf(f,'%7s  %8s  %9s\n',sta,lat,lon);
    end 
elseif Nc > 3 
    for iss = 1:NS
        sta = pad(stadir{iss,1},7);
        lat = pad(num2str(stadir{iss,2},'%7.4f'),8,'left');
        lon = pad(num2str(stadir{iss,3},'%8.4f'),9,'left');
        elv = pad(num2str(round(stadir{iss,4}),'%4d'),5,'left');
        fprintf(f,'%7s  %8s  %9s  %5s\n',sta,lat,lon,elv);
    end 
end
fclose(f);