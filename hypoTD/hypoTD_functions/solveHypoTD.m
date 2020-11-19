function [x,stats] = solveHypoTD(A,b,x0,max_iter,cgsteps,minZ,maxZ)
% function [x,stats] = solveHypoTD(A,b,x0,max_iter,cgsteps,minZ,maxZ)
%
% 2020-07-08
% Solves the hypoTD system using Newton's method. Allows depth bounds to
% be specified.
%
%   INPUTS
%
%        A = full hypoTD matrix
%        b = full hypoTD RHS vector
%       x0 = [E; N; Z] earthquake cartesian coordinates
% max_iter = maximum number of Newton iterations
%  cgsteps = number of Conjugate gradient steps to define Newton direction
%     minZ = minimum EQ depth (defaults to 0)
%     maxZ = maximum EQ depth (defaults to 1000)
%

if nargin < 7
    maxZ = 1000;
    if nargin < 6
        minZ = 0;
        if nargin < 5
            cgsteps = max(20,length(x)/2);
            if nargin < 4
                max_iter = 20;
            end
        end
    end
end


x  = x0;
L  = length(x);
Ne = L/3;
jz = 2*Ne+[1:Ne]';

msft  = norm(A*x-b);
CNVRG = 0;
iter  = 0;

max_try = 5;

while ~CNVRG && iter < max_iter
    
    dx = NewtonStepCG(A,b,x,cgsteps);
    
    for ii = 1:max_try
        
        xt = x + dx;
        xt(jz( xt(jz) < minZ )) = minZ;
        xt(jz( xt(jz) > maxZ )) = maxZ;
  
        msftT = norm(A*xt-b);

        if msftT < msft
            x = xt;
            msft = msftT;
            iter = iter + 1;
            break
        elseif ii < max_try
            dx = dx/2;
            continue
        else
            CNVRG = 1;
        end
    end
    
end
          
stats.iter = iter;
    
            
        
