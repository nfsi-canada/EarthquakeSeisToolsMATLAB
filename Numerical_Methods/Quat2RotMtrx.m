function R = Quat2RotMtrx(q)
% function R = Quat2RotMtrx(q)

% 2020-08-06
% Equation to convert a quaternion into a rotation matrix,
% taken from website:
% https://www.euclideanspace.com/maths/geometry/rotations/ ...
%             conversions/quaternionToMatrix/index.htm
%
% R = [1-2*(y^2+z^2)     2*(x*y-z*w)   2*(x*z+y*w)
%        2*(x*y+z*w)   1-2*(x^2+z^2)   2*(y*z-x*w)
%        2*(x*z-y*w)     2*(y*z+x*w) 1-2*(x^2+y^2)];
%
%  INPUTS 
%
%       q = [w,x,y,z] a quaternion
%
%  OUTPUTS
%
%       R = a 3x3 rotation matrix

R = [1-2*(q(3)^2+q(4)^2)           2*(q(2)*q(3)-q(4)*q(1))    2*(q(2)*q(4)+q(3)*q(1))
       2*(q(2)*q(3)+q(4)*q(1))   1-2*(q(2)^2+q(4)^2)          2*(q(3)*q(4)-q(2)*q(1))
       2*(q(2)*q(4)-q(3)*q(1))     2*(q(3)*q(4)+q(2)*q(1))  1-2*(q(2)^2+q(3)^2)];
       
