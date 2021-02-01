function G = buildStationDDmatrix(ev,g,v,W,Ne,swTD)
% function G = buildStationDDmatrix(ev,g,v,W,Ne,swTD)
%
% 2020-01-22
% This function builds a matrix for station-pari double-difference 
% hypocenter inversion
% 
% The matrix is structured [X Y Z (T)], i.e. the first Ne columns 
% correspond to the X-coordiante.
%
% INPUTS
%
%    evs == N x 1 vector of event indices 
%      g == N x 6 matrix of ray takeoff directions to the two stations
%           [x1 y1 z1, x2 y2 z2] 
%           for the two events (they should generally be very similar)
%      v == N x 1 vector of velocities at the source
%      W == N x 1 vector of weights for each equation
%     Ne == total number of events 
%   swTD == zero for double-difference system (default) 
%           or non-zero for triple-difference system
%
% OUTPUTS
%
%    G == N x (4 or 3)*Ne matrix with columns ordered [X Y Z (T)]
%         i.e. first ALL x-coordiantes, followed by ALL y

if nargin < 6
    swTD = 0;
end


N  = length(ev);
g1 = [g(:,1);  g(:,2);  g(:,3)];
g2 = [g(:,4);  g(:,5);  g(:,6)];
v  = repmat(v,3,1);

me = repmat(W,3,1).*(-(g1-g2)./v);  
j1 = repmat((1:N)',3,1);
j2 = [ev; ev+Ne; ev+2*Ne];

if swTD
    G = sparse(j1,j2,me,N,3*Ne);
else
    G = sparse(j1,j2,me,N,4*Ne);
end

    