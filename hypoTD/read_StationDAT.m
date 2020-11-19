function stadir = read_StationDAT(f)
% function stadir = read_StationDAT(f)
%
% 2020-10-09
% This function takes in a path to a hypoDD 'station.dat' formatted file
% and returns results in a cell array "stadir". If there is a 4th column
% it is treated as station elevation.
%
%   INPUTS
%
%      f == path to "station.dat"-style file
%
%   OUTPUTS
%
%     stadir == {STA,LAT,LON,(ELV)}
%
% 
A  = importdata(f);
ns = size(A.textdata,1);
nc = size(A.data,2);

stadir        = cell(ns,nc+1);
stadir(:,1)   = A.textdata;
stadir(:,2:nc+1) = num2cell(A.data(:,1:nc));
    
