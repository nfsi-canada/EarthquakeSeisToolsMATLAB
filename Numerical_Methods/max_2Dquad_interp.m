function [xmax,ymax,dmax] = max_2Dquad_interp(x,y,d)
%  function [xmax,ymax,dmax] = max_2Dquad_interp(x,y,d)
%
%  This function gets subsample precision of a maximum/minimum
%  by using a 2D quadratic interpolation around the sample-wise
%  maximum/mininum (and the eight adjacent points). It was intended 
%  for finding better cross-correlation delay times.
%  Notation for system is Ax = d. 
%
%  INPUT
%
%  x == the 3 x samples [dtmax-dt,dtmax,dtmax+dt]
%  y == the 3 y samples
%  d == 3x3 portion of the data 
%       (data can be e.g. the cross-correlation function)
%       corresponds to points [x-1,y-1
%                              x0 ,y-1
%                              x+1,y-1
%                              x-1,y0 ...]
%
%
%
%  OUTPUT 
%
%  xmax == corrected x-location of maximum
%  ymax == corrected y-location of maximum
%  dmax == estimated maximum cross-correlation

% -- Ensure input are column vectors
x  = mkcol(x);
y  = mkcol(y);
dv = d(:);

xv = repmat(x,3,1);
yv = repmat(y',3,1);
xv = xv(:);
yv = yv(:);

% -- Columns of A are (x^2, x, 1 -- multiplied by a,b,c)

% -- Equation begin modeled is   z = a*x^2 + b*y^2 + c*x + d*y + e*x*y + f

A = [xv.^2, yv.^2, xv, yv, xv.*yv, ones(9,1)];

% -- Solve quadratic for ax,ay,bx,by,c
abcdef = A\dv;

% -- Get x at maximum, new maximum of ccf;
a = abcdef(1);
b = abcdef(2);
c = abcdef(3);
d = abcdef(4);
e = abcdef(5);
f = abcdef(6);

den = 4*a*b-e^2;
if den > 0
    xmax = -(2*b*c-d*e)/den;
    ymax = -(2*a*d-c*e)/den;
    dmax = dot(abcdef,[xmax^2; ymax^2; xmax; ymax; xmax*ymax; 1]);
else
    xmax = NaN;
    ymax = NaN;
    dmax = NaN;
end

