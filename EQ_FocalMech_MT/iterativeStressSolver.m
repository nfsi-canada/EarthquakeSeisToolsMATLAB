function [x,msft] = iterativeStressSolver(A,b,tol,max_iter)
% function [x,msft] = iterativeStressSolver(A,b,tol,max_iter)
%
% 2020-08-07
% This is designed to account for a "flaw" in standard 
% linear stress inversions (using earthquake focal mechanisms).
% Because we don't know the amount of shear stress of each fault, we can't
% scale the RHS vector (b, composed of slip vectors) to scale, and assume
% unit length. This function iteratively solves and reweights the RHS 
% vector elements to account for this.
%
%  INPUTS
%
%     A = 3N x 5 matrix. Columns correspond to [T11 T12 T13 T22 T23]
%         First N rows correspond to east components, then north, then up.
%     b = 3N x 1 vector [east,north,up] components of observed slip vectors
%   tol = maximum difference between input and output slip vector lengths.
%         Default is 1e-3.
%   max_iter == Maximum number of iterations. Default is 50.
%
%  OUTPUTS
%
%     x = best-fit stress tensor elements [T11 T12 T13 T22 T23]'
%

if nargin < 4
    max_iter = 50;
    if nargin < 3
        tol = 1e-3;
    end
end

% -- The indices in 'jb' convert the RHS vector 'b' into a
% -- Nx3 matrix of slip vectors 'B'
N  = length(b)/3;
jb = [(1:N); N+(1:N); 2*N+(1:N)]';

% -- Initialize a vector of slip vector lengths
b0 = b;
nrmB = ones(N,1);

for ii = 1:50
    
    x = A\b;
    
    % -- Check the lengths of the predicted slip vectors
    bp = A*x;
    Bp = bp(jb);
    nrmBp = min(vecnorm(Bp')',1);
    
    if max(abs(nrmBp-nrmB)) > tol
        nrmB = nrmBp;
        b = b0.*repmat(nrmB,3,1);
    else
        break
    end
end

msft = norm(A*x-b)/norm(b);