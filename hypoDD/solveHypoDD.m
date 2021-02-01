function [x,stats] = solveHypoDD(A,b,z,max_iter,cgsteps,minZ,maxZ)
% function [x,stats] = solveHypoDD(A,b,z,max_iter,cgsteps,minZ,maxZ)
%
% 2022-01-12
% Solves the hypoDD system for optimal shifts using Newton's method. 
% Allows depth bounds to be specified.
%
%   INPUTS
%
%        A = full hypoDD matrix
%        b = full hypoDD RHS vector
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
%            [dE; dN; dZ; dT] earthquake cartesian coordinates
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
Ne = L/4;
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
        
        % -- Then correct times for projected events by moving all other
        % -- columns to RHS and solving
        jT = 3*Ne+[jz0; jz1];
        jC = setdiff(1:4*Ne,jT);
        bC = b-A(:,jC)*xt(jC);
        xt(jT) = A(:,jT)\bC;
        
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
stats.DT    = mean(abs(x(3*Ne+[1:Ne])));
stats.OS    = mean([x(1:Ne) x(Ne+[1:Ne]) x(2*Ne+[1:Ne]) x(3*Ne+[1:Ne])]);
            


% -- OLD VERSION WAS ACTUALLY UPDATING m INSTEAD OF SOLVING FOR dm!!!
%{

m  = m0;
L  = length(m);
Ne = L/4;
jz = 2*Ne+[1:Ne]';
jT = 2*Ne+[1:Ne]';
msft  = norm(G*m-d);
stats.msft0 = msft;
CNVRG = 0;
iter  = 0;

max_try = 5;

while ~CNVRG && iter < max_iter
    
    dm = NewtonStepCG(G,d,m,cgsteps);
    
    for ii = 1:max_try
        
        mt = m + dm;
        
        % -- Project depth if outside bounds    
        jz0 = find(mt(jz) < minZ);
        jz1 = find(mt(jz) > maxZ);
        mt(2*Ne+jz0) = minZ;
        mt(2*Ne+jz1) = maxZ;
        
        % -- Then correct times for projected events by moving all other
        % -- columns to RHS and solving
        jT = 3*Ne+[jz0; jz1];
        jC = setdiff(1:4*Ne,jT);
        dC = d-G(:,jC)*mt(jC);
        mt(jT) = G(:,jT)\dC;
        
        msftT = norm(G*mt-d);

        if msftT < msft
            m = mt;
            msft = msftT;
            iter = iter + 1;
            break
        elseif ii < max_try
            dm = dm/2;
            continue
        else
            CNVRG = 1;
        end
    end
    
end
          
stats.iter = iter;
stats.msft = msft;    
%}
        
