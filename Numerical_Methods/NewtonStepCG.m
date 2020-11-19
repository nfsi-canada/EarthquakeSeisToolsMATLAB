function [dx,flg] = NewtonStepCG(A,b,x,N)
% function [dx,flg] = NewtonStepCG(A,b,x,N)
%
% 2020-07-08
% This function computes a newton step for a A*x = b linear system, 
% but without ever computing the Hessian (A'*A), so that is stays fast 
% for large, sparse systems.
%
%  INPUTS
%
%     A = matrix (A*x = b)
%     b = RHS vector 
%     x = initial estimate of x
%     N = maximum number of steps for conjugate-gradient method
%
%  OUTPUT
%
%    dx = Newton step (x1 = x0 + dx)
%   flg = 0 if things went well, 1 if there's an issue

% CG Algorithm comes from Nocedal and Wright textbook 
% page 112, algorithm 5.2 ...  their notation = my notation:
% A = A'*A
% b = g
% x = dx
% p = ddx

% -- Steepest descent direction
ATb = A'*b;
g = -(A'*(A*x)-ATb);
L = length(x);

% -- Initialize CG descent direction (dx)
dx  = zeros(length(x),1);
r   = A'*(A*dx)-g;
ddx = -r;

% -- Conjugate-gradient (CG) method for step direction
for ii = 1:min(N,L);
    
    ATAddx = A'*(A*ddx);
    
    a  = (r'*r)/(ddx'*ATAddx); % Alpha - a step length
    dx = dx + a*ddx;           % Actually changing dx here     
    ru  = r + a*ATAddx;       % Find new residual
    bb  = (ru'*ru)/(r'*r);    % Beta - a weight of new vs old direction
    ddx = -ru + bb*ddx;       % Update step direction
    r   = ru;                 % Update 'old' residual term
    
end

% -- Step-length
ATAdx = A'*(A*dx);
f1    = x'*ATAdx - dx'*ATb;
f2    = dx'*ATAdx;
dx    = -(f1/f2)*dx;

if f1 > 0 
    flg = 1;    
else
    flg = 0;
end



