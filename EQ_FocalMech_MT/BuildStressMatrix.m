function A = BuildStressMatrix(n)
% function A = BuildStressMatrix(n)
%
% 2020-08-06
% This function builds a matrix to invert for the state of tectonic
% stress using Vavrycuk's [2014,GJI,Eq 8] notation (but not transposed).
% 
% For a compression=positive sign convention
%  
%   INPUTS
%
%    n = N x 3 matrix of fault normals [n1,n2,n3]
%        n3 MUST BE POSITIVE (or zero)!! AND THIS DOESN'T CHECK THAT.
%        (normal must point up)
%
%   OUTPUTS
%
%    A = (Nx3) x 5 matrix to invert for stress
%        Row order is all East/s1 equations, then North/s2, then Up/s3
%        Column order is [S11, S12, S13, S22, S23]

n1 = n(:,1);
n2 = n(:,2);
n3 = n(:,3);

A = [n1.*(n2.^2+2*n3.^2), n2.*(1-2*n1.^2), n3.*(1-2*n1.^2), n1.*(-n2.^2+n3.^2),  -2*n1.*n2.*n3,
     n2.*(-n1.^2+n3.^2),  n1.*(1-2*n2.^2),-2*n1.*n2.*n3,    n2.*(n1.^2+2*n3.^2),  n3.*(1-2*n2.^2),
     n3.*(-2*n1.^2-n2.^2),-2*n1.*n2.*n3,   n1.*(1-2*n3.^2), n3.*(-n1.^2-2*n2.^2), n2.*(1-2*n3.^2)]; 
 