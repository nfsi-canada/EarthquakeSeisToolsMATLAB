function x1 = MovingMean(x0,N)
% function x1 = MovingMean(x0,N)
%
%
% 2020-03-28
% This function produces a moving average over N values for each column in the 
% input signal matrix x0

[N1,N2] = size(x0);

% -- Build matrix A to perform averaging
% -- Slightly different if N is even than is N is odd
% -- First get indices of points to include in each moving average j1 to j2
A = zeros(N1);
if N/2-floor(N/2) > 0.2
    h  = fix((N-1)/2);
    j1 = [1:N1]'-h;
    j2 = [1:N1]'+h;
else
    h  = fix(N/2);
    j1 = [1:N1]-h;
    j2 = [1:N1]+h-1;
end
j1(j1<1)  = 1;
j2(j2>N1) = N1;

% -- Now actually build averaging matrix row by row
for ii = 1:N1
    A(ii,j1(ii):j2(ii)) = 1/(j2(ii)-j1(ii)+1);
end

% -- Perform averaging with matrix multiplication
x1 = A*x0;



