function [m]=sdr2mij(s,d,r)
% function [m]=sdr2mij(s,d,r)
%
% Converts strike S, dip D, rake R (in degrees) 
% to elements of moment tensor using Aki & Richards
% formulae followed by transformation from N,E,D to 
% E,N,U.

% Move to radians.
d2r=pi/180;
ss=s*d2r;
dd=d*d2r;
rr=r*d2r;

% Aki and Richard (N,E,D) convention.
m=zeros(3,3);
m(1,1)=-(sin(dd)*cos(rr)*sin(2*ss)+sin(2*dd)*sin(rr)*sin(ss)^2);
m(1,2)=(sin(dd)*cos(rr)*cos(2*ss)+0.5*sin(2*dd)*sin(rr)*sin(2*ss));
m(1,3)=-(cos(dd)*cos(rr)*cos(ss)+cos(2*dd)*sin(rr)*sin(ss));
m(2,2)=(sin(dd)*cos(rr)*sin(2*ss)-sin(2*dd)*sin(rr)*cos(ss)^2);
m(2,3)=-(cos(dd)*cos(rr)*sin(ss)-cos(2*dd)*sin(rr)*cos(ss));
m(3,3)=sin(2*dd)*sin(rr);
m(2,1)=m(1,2);
m(3,1)=m(1,3);
m(3,2)=m(2,3);

% Rotate into SW3D model convention (E,N,U)
rotm=[0 1 0; 1 0 0; 0 0 -1];
m=rotm*m*rotm';

return


