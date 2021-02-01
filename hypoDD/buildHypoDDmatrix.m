function G = buildHypoDDmatrix(evs,g,v,W,Ne)
% function G = buildHypoDDmatrix(evs,g,v,W,Ne)
%
% 2020-01-07
% This function builds a matrix for double-difference hypocenter inversion
% following equations 6-9 from Waldhauser and Ellsworth [2000,BSSA].
% I intended this to build the matrix for either CT or CC data at once, so 
% you may need to call it twice and combine the two G.
% 
% The matrix is structured [X Y Z T], i.e. the first Ne columns correspond
% to the X-coordiante.
%
% INPUTS
%
%    evs == N x 2 lists of event numbers 
%      g == N x 6 matrix of [x1 y1 z1, x2 y2 z2] ray takeoff directions 
%           for the two events (they should generally be very similar)
%      v == N x 2 matrix [v1, v2] of velocities at the sources
%      W == N x 1 vector of weights for each equation
%     Ne == total number of events 
%
% OUTPUTS
%
%    G == N x 4*Ne matrix with columns ordered [X Y Z T]
%         i.e. first ALL x-coordiantes, followed by ALL y

evA = evs(:,1);
evB = evs(:,2);
N   = length(evA);

gA = [g(:,1); g(:,2); g(:,3)];
gB = [g(:,4); g(:,5); g(:,6)];
vA = repmat(v(:,1),3,1);
vB = repmat(v(:,2),3,1);


% -- Should the spatial parts be divided by 2??!!!??
me =  repmat(W,8,1).*[-gA./vA; ones(N,1); gB./vB; -ones(N,1)];  
j1 = repmat([1:N]',8,1);
j2 = [evA; evA+Ne; evA+2*Ne; evA+3*Ne; evB; evB+Ne; evB+2*Ne; evB+3*Ne];

% -- Add zero-shift constraint if W0 is given
G = sparse(j1,j2,me,N,4*Ne);

    