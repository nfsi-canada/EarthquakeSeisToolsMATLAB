function [FM,SNR] = FirstMotionSNR(x,p,n);
% function [FM,SNR] = FirstMotionSNR(x,n);
%
% 2020-01-12
% This function gets the first-motion of a P-wave and computes a 
% signal-to-noise ratio (SNR). 
%
%   INPUTS
%
%        x == the P-waveform, should include at least 'n' oscillations
%             prior to P
%        p == index of the P-pick in x
%        n == number or pre-P peaks to take maximum of for SNR
%
%   OUTPUTS
%
%       FM == the first-motion, +1/-1 for up/down
%      SNR == amplitude of first P-pulse divided by avearage
%             of 'n' previous peaks

L = length(x);
x = x/max(abs(x)); % -- Not necessary unless plotting!!

% -- Zero-crossings of waveform
s = sign(x);
d = [0; diff(s)];
z = [1; find(d); L];

% -- Check if P-pick is before or after "it's peak", meaning
% -- the (absolute) maxima between zero-crossings
% -- or what if it is the peak???

j0 = find(z <= p,1,'last');
j1 = find(z  > p,1);
z0 = z(j0);
z1 = z(j1);
[ampP,jm] = max(abs( x(z0:z1-1) ));
mx = z0+jm-1;

if mx <= p
    j0 = j1;
    z0 = z1;
    j1 = j0+1;
    z1 = z(j1);
    ampP = max(abs( x(z0:z1-1) ));
end

% -- Get first-motion
FM = sign(x(z0));

% -- For pre-P peaks
ampN = zeros(n,1);
for ii = 1:min(n,j0-1)
    j1 = j0;
    z1 = z0;
    j0 = j1-1;
    z0 = z(j0);
    ampN(ii) = max(abs( x(z0:z1-1) ));
end

% -- Get signal-to-noise ratio
SNR = ampP/max(ampN);

%{
figure(1)
clf 
hold on
plot([1:L],x/max(abs(x)))
plot([p,p],[-1,1],'r--')
xlim([p-50,p+50])
grid on
keyboard
%}