function G = buildHypoTDmatrix(evs,g,v,W,Ne)
% function G = buildHypoTDmatrix(evs,g,v,W,Ne)
%
% 2020-01-15
% This was modified from "buildHypoDDmatrix" to build a matrix for
% triple-difference hypocenter inversion.
%
% 2020-01-07
% This function builds a matrix for double-difference hypocenter inversion
% following equations 6-9 from Waldhauser and Ellsworth [2000,BSSA].
% I intended this to build the matrix for either CT or CC data at once, so 
% you may need to call it twice and combine the two G.
% 
% The matrix is structured [X Y Z], i.e. the first Ne columns correspond
% to the X-coordiante.
%
% INPUTS
%
%    evs == N x 2 lists of event numbers 
%      g == N x 12 matrix of ray takeoff directions 
%           [xA1 yA1 zA1, xB1 yB1 zB1, xA2 yA2 zA2, xB2 yB2 zB2] 
%           for the two events (they should generally be very similar)
%      v == N x 2 matrix [vA, vB] of velocities at the sources
%      W == N x 1 vector of weights for each equation
%     Ne == total number of events 
%
% OUTPUTS
%
%    G == N x 3*Ne matrix with columns ordered [X Y Z]
%         i.e. first ALL x-coordiantes, followed by ALL y

evA = evs(:,1);
evB = evs(:,2);
N   = length(evA);

gA1 = [ g(:,1);  g(:,2);  g(:,3)];
gB1 = [ g(:,4);  g(:,5);  g(:,6)];
gA2 = [ g(:,7);  g(:,8);  g(:,9)];
gB2 = [g(:,10); g(:,11); g(:,12)];
vA = repmat(v(:,1),3,1);
vB = repmat(v(:,2),3,1);


me = repmat(W,6,1).*[-(gA1-gA2)./vA; (gB1-gB2)./vB];  
j1 = repmat([1:N]',6,1);
j2 = [evA; evA+Ne; evA+2*Ne; evB; evB+Ne; evB+2*Ne];
G  = sparse(j1,j2,me,N,3*Ne);

    