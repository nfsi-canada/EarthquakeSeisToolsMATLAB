function V = Vfun(X,FV,Vmin,Vmax)
% function V = Vfun(X,FV,Vmin,Vmax)
% 
% 2020-12-18
% Interpolates velocity model FV at 'X' and ensures output is between
% Vmin and Vmax. Mainly for use in RayTrace3D.
%
%  INPUTS
%
%     X == [N x 3] array of points to sample FV at x,y,z 
%    FV == interpolant of velocity, made with 'scatteredInterpolant'
%          or 'griddedInterpolant'. It takes in either [x,y,z] locations.
%  Vmin == project interpolated velocities < Vmin to Vmin
%  Vmax == project interpolated velocities > Vmax to Vmax
%
%  OUTPUTS
%
%    V == velocities at points "X"


V = FV(X);
V(V<Vmin) = Vmin;
V(V>Vmax) = Vmax;
