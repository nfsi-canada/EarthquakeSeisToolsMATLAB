function [x,stats] = solveHypoTD(A,b,z,max_iter,cgsteps,minZ,maxZ)
% function [x,stats] = solveHypoTD(A,b,z,max_iter,cgsteps,minZ,maxZ)
%
% 2022-01-12
% Solves the hypoTD system for optimal shifts using Newton's method. 
% Allows depth bounds to be specified.
%
% This was modified from "solveHypoDD." to work without columns for time.
%
%   INPUTS
%
%        A = full hypoTD matrix
%        b = full hypoTD RHS vector
%        z = initial solution depths (only needed for depth constraints)
% max_iter = maximum number of Newton iterations
%  cgsteps = number of Conjugate gradient steps to define Newton direction
%     minZ = minimum EQ depth (defaults to 0)
%     maxZ = maximum EQ depth (defaults to 1000)
%
% OUTPUTS 
%
%       x == optimal shift in hypocenters "dm", but respects earthquake
%            depth bounds
%            [dE; dN; dZ] earthquake cartesian coordinates
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

L  = size(A,2);
Ne = L/3;
x  = zeros(L,1);
jz = 2*Ne+(1:Ne)';

msft  = norm(A*x-b);
stats.msft0 = msft;
CNVRG = 0;
iter  = 0;
max_try = 5;

while ~CNVRG && iter < max_iter
    
    dx = NewtonStepCG(A,b,x,cgsteps);
    
    for ii = 1:max_try
        
        xt = x + dx;
        
        % -- Project depth if outside bounds 
        jz0 = find(z+xt(jz) < minZ);
        jz1 = find(z+xt(jz) > maxZ);
        xt(2*Ne+jz0) = minZ-z(jz0);
        xt(2*Ne+jz1) = maxZ-z(jz1);
                
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
          
stats.iter  = iter;
stats.msft  = msft;    
stats.prctM = 100*(stats.msft-stats.msft0)/stats.msft0;
stats.DH    = mean(sqrt(x(1:Ne).^2 + x(Ne+[1:Ne]).^2));
stats.DZ    = mean(abs(x(2*Ne+[1:Ne])));
stats.OS    = mean([x(1:Ne) x(Ne+[1:Ne]) x(2*Ne+[1:Ne])]);
            
