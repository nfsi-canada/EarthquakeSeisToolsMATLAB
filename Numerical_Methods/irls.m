function [m] = irls(A,d,p)

% FUNCTION M = IRLS(A,D,P)
% Function to compute iteratively reweighted least squares
% solution to A m = d, for a p norm.
    
% Initialization.
A0=zeros(size(A));
d0=zeros(size(d));

% First iteration solve using least squares.
m=A\d;
%     return

% While loop with convergence criterion.
eps=0.0001;
phi=1000.0;
k=0;
while phi > 0.005 
    m0=m;
    k=k+1;
    disp(['Iteration: ',num2str(k)])

    % Compute residual vector, and take those near zero values
    % to epsilon..
    res=abs(A*m0-d);
    ix=find(res < eps);
    res(ix)=eps;
    res=res.^(p-2);
    res=res/max(res);

    % Now convert to diagonal weighting elements according to p-norm.
    % Extra square root is just to allow the problem to be solved using
    % the matlab (least squares) division, ie A'*R*(A*m-d)=0 or
    % A'*sqrt(R)'*(sqrt(R)*A*m-sqrt(R)*d)=B'*(B*m-d0)=0. Thus just
    % solve a modified system where A-->B=sqrt(R)*A and d-->d0=sqrt(R)*d.
    res=sqrt(res);

    % Create equivalent least-squares problem.
    d0=res.*d; 
    for i=1:size(A,1);
        A0(i,:)=res(i)*A(i,:);
    end
    m=A0\d0;
    phi=sqrt((m0-m)'*(m0-m));
end
