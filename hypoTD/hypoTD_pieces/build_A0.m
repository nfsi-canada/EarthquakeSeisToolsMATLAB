% 2020-05-06
% This piece of hypoTD builds A0 and b0, which set the solution equal to 
% the current solution. This part of the system is meant to be weighted 
% low, but high enough to enusre no hypocenter coordinates are 
% underdetermined. 
%
% After the inversion, the matrix can be analyzed to see how important this
% constraint was to each event.

A0 = speye(3*Ne);
b0 = x0;
.