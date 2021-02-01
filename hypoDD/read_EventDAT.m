function evdir = read_EventDAT(f)
% function evdir = read_EventDAT(f)
%
% 2020-10-09
% This function takes in a path to a hypoDD 'event.dat' formatted file
% and returns results in a cell array "evdir"
%
%   INPUTS
%
%      f == path to "event.dat"-style file
%
%   OUTPUTS
%
%     evdir == {DATETIME,LAT,LON,DEP,MAG}
%
A     = importdata(f);
ne    = size(A,1);
evdir = cell(ne,6);
for ie = 1:ne
    date = num2str(A(ie,1));
    time = num2str(A(ie,2)+1e8);
    evdir{ie,6} = strcat(date,time(2:end-2));
    evdir{ie,1} = datetime(strcat(date,time(2:end)),'InputFormat','yyyyMMddHHmmssSS');
end
evdir(:,2:5) = num2cell(A(:,3:6));