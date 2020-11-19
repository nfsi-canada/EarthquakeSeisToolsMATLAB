function x = gen_resample(x,newdt,olddt);
% function x = gen_resample(x,newdt,olddt);
%
% 2020-07-17
% This function takes in a signal, an old sample interval,
% and a new sample interval. It then uses MATLAB's built-in "resample"
% function, after using "rat" computing the integer input values
% needed for resample.
%
%  INPUTS
%
%       x = original input signal
%   newdt = desired sample interval (generally in seconds)
%   olddt = desired sample interval (units matching newdt)
%
%  OUTPUTS
%
%     x = resampled version of original signal

r = olddt/newdt;
[p,q] = rat(r);
x = resample(x,p,q);

