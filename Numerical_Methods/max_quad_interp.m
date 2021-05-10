function [xmax,dmax] = max_quad_interp(x,d)
%  function [xmax,dmax] = max_quad_interp(x,d)
%
%  This function gets subsample precision of a maximum/minimum
%  by using a quadratic interpolation around the sample-wise
%  maximum/mininum (and the two adjacent points). It was intended 
%  for finding better cross-correlation delay times.
%  Notation for linear system is Ax = d. 
%
%  2019-02-22
%  Added check to ensure interpolation, no extrapolation...
%
%  INPUT
%
%  x == the 3-sample time vector [dtmax-dt,dtmax,dtmax+dt]
%  d == 3 sample long portion of the data
%       (data can be e.g. the cross-correlation function)
%
%  OUTPUT 
%
%  xmax == corrected maximum time
%  dmax == estimated maximum cross-correlation

% -- Ensure input are column vectors
x = x(:);
d = d(:);

% -- Columns of A are (x^2, x, 1 -- multiplied by a,b,c)
% -- for a  ax^2 + bx + c = 0 model
A = [x.^2, x, ones(3,1)];

% -- Solve quadratic for a,b,c
abc = A\d;

% -- Get x at maximum, new maximum of ccf;
% -- Ensure the maxima (or minima) is within these three samples
xmax = -abc(2)/(2*abc(1));
xmax = min( max(xmax,min(x)),max(x) );

% -- Get max (or min) of d
dmax = dot(abc,[xmax.^2; xmax; 1]);
    
    
    
