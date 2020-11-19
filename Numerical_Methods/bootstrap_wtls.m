function [slp,b,bnd] = bootstrap_wtls(N,ci,x,y,ex,ey)
% function [slp,b,bnd] = bootstrap_wtls(N,ci,x,y,ex,ey)
%
% This function does bootstrap method to get error estimates on the
% slope for a total-least squares line regression.   
%
%     INPUTS
%
%         N == number of times to run experiment
%        ci ==  confidence interval (multiple of S.D.) -- must be 1,2, or 3
%      x, y == data to to be fit
%     ex,ey == individual error esitmates for each point, optional
%
%     OUTPUTS
%
%     slp == mean slope %the resulting slope using all the data
%       b == mean y-intercept %the y-intercept of the regression using all the data 
%     bnd == upper and lower bound of slope based of ci
%

% -- Get number of points
np = length(x);

% -- First check if errors were provided
if nargin < 5
	ex = ones(np,1);
	ey = ones(np,1);
elseif length(ex)==1
	ex = ones(np,1);
	ey = ones(np,1);
end 


% -- Create array to store all the slopes. 
a = zeros(N,1);
b = zeros(N,1);


% -- Do initial regression with full dataset
%[slp,b] = robust_wtls(x,y,ex,ey,3,50);

for ii = 1:N

	inds = randsample(np,np,true);
	rx  = x(inds);
	ry  = y(inds);
	rex = ex(inds);
	rey = ey(inds);

	% -- Use robust weighted total-least-squares algorithm
	% -- Eliminate points < 3 s.d. away from line, use 50 iterations
	[a(ii),b(ii)] = robust_wtls(rx,ry,rex,rey,3,50);

end

% -- Take mean of distributions...
slp = mean(a);
b   = mean(b);

% -- Calculate error at desired confidence interval
% -- Start by selecting percentiles based of confidence interval
if ci == 1
	lp = 15.9;
	up = 84.1;
elseif ci == 2
	lp = 2.3;
	up = 97.7;
elseif ci == 3
	lp = 0.1
	up = 99.9
else
	disp('Error: C.I. not 1,2, or 3 in bootstrap')
end


% -- Sort 'a' and interpolate for bounds
a   = sort(a); 
p   = [0:100/(N-1):100];
bnd = interp1(p,a,[lp,up]);




% -- Get standard deviation, then bounds
%sa  = std(a);
%bnd = [slp-ci*sa; slp+ci*sa];
%err = ci*std(a);

