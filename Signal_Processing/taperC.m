function [X] = taperC(X,dt,dur1,dur2,t1,t2);
% function [X] = taperC(X,dt,dur1,dur2,t1,t2);
%  
% 2019-02-01
% This function was modified from MGB's code taper. 
%
%  INPUTS
%
%         X == one or more column vectors
%        dt == the sample interval (in seconds)
%      dur1 == the duration of the first (ramp up) taper, in seconds
%      dur2 == duration of ramp-down taper (if not provided, dur1 is used)
%        t1 == time to end ramp-up, in sec. If not provided, dur1 is used.
%              Must be >= dur1.
%        t2 == time to begin ramp-down, in sec. If not provided, 
%              T-dur2 is used. Must be <= T-dur2;
%
%  OUTPUTS
%
%        X == tapered waveforms

[nx,nc] = size(X);
taper   = ones(nx,1);

if nargin < 6;   
    if nargin < 5
        t1 = dur1;
    end
    if nargin < 4
        dur2 = dur1;
    end
    
    t2 = nx*dt-dur2;
end

if t1 > 0
    N  = fix(dur1/dt);
    ii = [0:N]'*dt/dur1;
    ct = 0.5-0.5*cos(pi*ii);
    jj = fix(t1/dt+1);
    
    taper(jj-N:jj) = ct;
    taper(1:jj-N)  = 0;
end
if (t2 > 0)*(t2 <= nx*dt-dur2)
    
    N  = fix(dur2/dt);
    ii = [0:N]'*dt/dur2;
    ct = 0.5+0.5*cos(pi*ii);
    jj = fix(t2/dt);
    
    taper(jj:jj+N) = ct;
    taper(jj+N:nx) = 0;
end

% -- Apply taper to each column of X
X = repmat(taper,1,nc).*X;


