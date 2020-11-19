function q = RandomQuaternion(N)
% function q = RandomQuaternion(N)
%
% 2020-08-06
% This function generates (uniform) random quaternions (w,x,y,z) with 
% the algorithm from Kuffner 2004,IEEE Conf. on Robotics and Automation
% A single quaternion can be converted to a rotation matrix with the 
% function Quat2RotMtrx
%
% INPUTS
%
%    N == simply an integer number of quaternions to produce
%
% OUTPUTS
%
%    q = [w,x,y,z] an N x 4 matrix of random (unit) quaternions

q = zeros(N,4);
s = rand(N,1);
s1 = sqrt(1-s);
s2 = sqrt(s);
t1 = 2*pi*rand(N,1);
t2 = 2*pi*rand(N,1);
q(:,1) = cos(t2).*s2;
q(:,2) = sin(t1).*s1;
q(:,3) = cos(t1).*s1;
q(:,4) = sin(t2).*s2;