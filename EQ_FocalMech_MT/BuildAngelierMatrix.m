function A = BuildAngelierMatrix(n,v)
% function A = BuildAngelierMatrix(n,v)
%
% 2020-08-06
% This function builds a matrix to invert for the state of tectonic
% stress using a method similar to Angelier [2002,GJI]
%
% The point is we are trying to maximize the dot(T*n,v)==dot(T*v,n)...
% 
% 
%  
%   INPUTS
%
%    n = N x 3 matrix of fault normals [n1,n2,n3]
%        n3 MUST BE POSITIVE (or zero)!! AND THIS DOESN'T CHECK THAT.
%        (normal must point up)
%
%   OUTPUTS
%
%    A = N x 5 matrix to invert for stress
%        Column order is [S11, S12, S13, S22, S23]

n1 = n(:,1);
n2 = n(:,2);
n3 = n(:,3);
v1 = -v(:,1);
v2 = -v(:,2);
v3 = -v(:,3);

A = sum([n1.*v1-n3.*v3, n1.*v2+n2.*v1, n1.*v3+n3.*v1, n2.*v2-n3.*v3, n2.*v3+n3.*v2]);

          
% -- What would the RHS equal... trying to maximize it...


% -- Maybe I can use a descent method to keep the stress tensor
% -- have a "scalar moment" of 1?
%M0 = sqrt(S11^2 + S12^2 + S13^2 + S22^2 + S23^2 + S11*S22 )
%M0 = (1/sqrt(2))* sqrt(S11^2 + 2*S12^2 + 2*S13^2 + S22^2 + 2*S23^2 + (-S11-S22)^2);
% rM0 = sqrt(S11^2 + 2*S12^2 + 2*S13^2 + S22^2 + 2*S23^2 + S11^2 + S22^2 + 2*S11*S22 );
% Ignore square root!
%{
M0 = S11^2 + S12^2 + S13^2 + S22^2 + S23^2 + S11*S22;
JM = zeros(5,1);
JM(1) = 2*S11+S22;
JM(2) = 2*S12;
JM(3) = 2*S13;
JM(4) = 2*S22+S11;
JM(5) = 2*S23;

HM = diag(2*ones(5,1));
HM(1,4) = 1;
HM(4,1) = 1;
%}
