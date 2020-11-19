function M = sdrs2mij(sdrs,vpvs)
% function M = sdrs2mij(sdrs,vpvs)
%
% This function calculates a moment tensor from a strike, dip, rake, and slope
% (s,d,r,p). This follows Vavrycuk 2011 JGR. I assume the source region is 
% isotropic and a Poisson solid unless a vpvs or Poisson ratio is provided.
%
% Vavrycuk has an equation that can get lambda,mu ratio (and thus vp/vs) 
% from M1,M2,M3...
%
% NEED TO FIGURE OUT WHICH CONVENTION TO USE.
%
%    INPUTS
%
%       sdrs = [s,d,r,a]'
%       s = fault strike usindg right-hand-rule (0-360 clockwise from N)
%       d = fault dip usindg right ahnd rule (0-90 from horizontal)
%       r = fault rake, describes sense of motion, describes motion of upper
%           block relative to strike. In degrees. 
%               0       = left-lateral
%               +/- 180 = right-lateral
%               +90     = reverse
%               -90     = normal
%       a = slope, in degrees, angle from fault plane of dislocation vector
%           0 = shear, +90 = extensional, -90 = compresional
%
%    vpvs = Vp/Vs if > 1, Poisson's ratio if < 0.5, if not provided, default is
%           Vp/Vs = 1.73 (Poisson's ratio =  0.25)
%
%    OUTPUTS
%
%      M = the normalized moment tensor with E,N,U orientation. 

% -- Calculate mu/lamda ratio mlr
if nargin < 2
    vpvs = sqrt(3);
elseif vpvs < 0.5
    vpvs = poisson2vpvs(vpvs);
end
mlr = 1/(vpvs^2-2);

% -- Decompose sdrs
s = sdrs(1);
d = sdrs(2);
r = sdrs(3);
a = sdrs(4);


% -- Calculate fault normal and dislocation vector
n1 = -sind(d)*sind(s);
n2 =  sind(d)*cosd(s);
n3 = -cosd(d);

v1 = (cosd(r)*cosd(s) + cosd(d)*sind(r)*sind(s))*cosd(a) - sind(d)*sind(s)*sind(a);
v2 = (cosd(r)*sind(s) - cosd(d)*sind(r)*cosd(s))*cosd(a) + sind(d)*cosd(s)*sind(a);
v3 = -sind(r)*sind(d)*cosd(a) - cosd(d)*sind(a);

% -- Potency tensor D without scaling factors
D = [ 2*n1*v1    ,  n1*v2+n2*v1,  n1*v3+n3*v1
      n1*v2+n2*v1,    2*n2*v2  ,  n2*v3+n3*v2
      n1*v3+n3*v1,  n2*v3+n3*v2,      2*n3*v3];

% -- Form M assuming Isotropic, Poisson Solid (lambda = mu)
% -- If I wanted to allow different poisson's ratio, note that:
% -- mlr == mu/lambda = ((vp/vs)^2 - 2)^-1
M = eye(3)*trace(D) + 2*mlr*D;

% -- Normalize M to a moment of 1 N m 
Mo = sqrt(sum(M(:).^2))/sqrt(2);
M  = M/Mo;

% -- Rotate from N,E,D to E,N,U?
rm = [0 1 0; 1 0 0; 0 0 -1];
M  = rm*M*rm';





% Aki and Richard (N,E,D) convention. FOR DOUBLE-COUPLE!!!!
%M=zeros(3,3);
%M(1,1) = -(sind(d)*cosd(r)*sind(2*s) + sind(2*d)*sind(r)*sind(s)^2);
%M(1,2) =  (sind(d)*cosd(r)*cosd(2*s) + 0.5*sind(2*d)*sind(r)*sind(2*s));
%M(1,3) = -(cosd(d)*cosd(r)*cosd(s)   + cosd(2*d)*sind(r)*sind(s));
%M(2,2) =  (sind(d)*cosd(r)*sind(2*s) - sind(2*d)*sind(r)*cosd(s)^2);
%M(2,3) = -(cosd(d)*cosd(r)*sind(s)   - cosd(2*d)*sind(r)*cosd(s));
%M(3,3) =  (sind(2*d)*sind(r));
%M(2,1) = M(1,2);
%M(3,1) = M(1,3);
%M(3,2) = M(2,3);
