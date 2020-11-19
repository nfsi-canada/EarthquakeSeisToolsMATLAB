function [bp,f,P1,P0] = OptBand(x1,x0,dt,fmin,fmax)
% function [bp,f,P1,P0] = OptBand(x1,x0,dt,fmin,fmax)
%
% 2019-02-01
% This function takes in a target waveform and a background waveform,  
% and computes an optimal bandpass (to maximize SNR in the target waveform). 
%
% For applications that compare waveforms (e.g. alignment), you probably should
% not use frequencies as high as HF, because although it has high SNR, the
% Green's functions for two events may not be similar at such high frequencies.
% In most cases, HF = fmax. 
%
%
%     INPUTS
%
%         x1 = window with the target signal 
%         x0 = window for comparison (e.g. pre-signal background noise)
%              Note: if including multiple components they should be provided 
%                    in separate columns of x0,x1
%         dt = the sample interval, in seconds. 
%       fmin = minimum frequency of relevance (default 0.1 Hz)
%       fmax = maximum frequency of relevance (default 15 Hz)
%
%
%    OUTPUTS
% 
%         hp = suggested highpass (in Hz)         
%          f = the corresponding frequencies for P0,P1 (logarithmic spacing). 
%         P1 = power spectrum of target signal
%         P0 = background power spectrum
%       %%%SR = the ratio of normalized power spectra 
%

if nargin < 5
    fmax = 15;
    if nargin < 4
        fmin = 0.1;
    end
end

L1 = length(x1);
L0 = length(x0);
T1 = L1*dt;
T0 = L0*dt;
LM = max(L1,L0);

% -- First taper both waveforms
x1 = taperC(x1,dt,T1/12);
x0 = taperC(x0,dt,T0/12);

% -- Normalize power in each waveform
%x1 = x1/norm(x1);
%x0 = x0/norm(x0);

% -- Pad each to the same length (next power of two). 
% -- Can just do in fft command instead of beforehand
LC = 2^ceil(log(LM)/log(2)); 
x1 = pad0(x1,LC);
x0 = pad0(x0,LC);

% -- Column-vector of corresponding frequencies (in Hz)
fe  = [0:LC/2]'/(LC*dt);

% -- Resample to frequencies logarithmically...
lgf = linspace(log10(fe(2)),log10(fe(end)),LC/2)';
f   = 10.^lgf;

% -- Compute power spectra using multitaper estimate
% -- Sum across components (E,N,Z)?
% -- Correct relative scale by for length of signals included
P1 = sum(pmtm(x1),2);
P0 = sum(pmtm(x0),2);
P1 = P1*L0/L1;

% -- Sample spectra at log-spaced frequencies
P1 = interp1(fe,P1,f,'pchip','extrap');
P0 = interp1(fe,P0,f,'pchip','extrap');

% -- I want to select a continous band of frequencies with the largest 
% -- difference in the signal energy - noise energy
[~,jmin] = min(abs(f-fmin));
[~,jmax] = min(abs(f-fmax));
N        = jmax-jmin+1;

% -- Finding optimal band, intergrate P1-P0 in log space, from each point
% -- Subtract 0.3 from dP (~factor of 2)
dP = log10(P1(jmin:jmax))-log10(P0(jmin:jmax))-0.3;
C  = flipud(cumsum(flipud(dP)));
 
[mx,jb] = max(C);

% -- Because sometimes this was selecting quite low frequencies, do another 
% -- check on low end. This avoids including a wide band of low-frequencies
% -- with SNR barely above 1. 
mdP = mean(dP(jb:jmax));
jb  = jb+find(dP(jb:jmax)>0.5*mdP,1)-1;

% -- Frequencies of band with optimal SNR
hp = f(jb+jmin-1);
hp(hp < fmin) = fmin;
hp(hp > fmax) = fmax;

disp('THIS HAS NOT BEEN TESTED. 2020-08-20')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Plotting
%{
% -- Get log spectral ratio
SR = log10(P1./P0);

% -- What if we take only positive log-ratio values, fit parabola.  
jp  = find(SR);
A   = [f(jp).^2, f(jp), ones(length(jp),1)];
d   = SR(jp);
abc = A\d;
rSR = abc(1)*f.^2 + abc(2)*f + abc(3);
semilogx(f,rSR,'r');

sfigure(444);
clf
subplot(211)
loglog(f,P0,'k',f,P1,'r')
hold on
loglog([bp(1),bp(1)],[min(P0),max(P0)],'g--',[bp(2),bp(2)],[min(P0),max(P0)],'c--')
grid on
ylabel('Power Spectra')

subplot(212)
semilogx(f,SR,'k',f,rSR,'r')
hold on
plot([bp(1),bp(1)],[min(SR),max(SR)],'g--',[bp(2),bp(2)],[min(SR),max(SR)],'c--')
grid on
xlabel('Frequency (Hz)')
ylabel('Log Power Ratio')

fprintf('Opt. BP %4.1f to %4.1f Hz\n',bp(1),bp(2));
keyboard
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Old stuff

% -- Normalize them to the same total power
%P1 = P1/sum(P1);
%P0 = P0/sum(P0);

% -- Looking at the semilogx plot (rSR vs. f), find points of maximum upward 
% -- curvature, and the inflection point to the right.  
%drdf = diff(rSR);
%d2df = diff(drdf);

%[~,j2] = max(drdf);
%[~,j1] = max(d2df);


%f([j1,j2])
%plot([f(j1),f(j1)],[min(SR),max(SR)],'g--',[f(j2),f(j2)],[min(SR),max(SR)],'c--')
% -- How about I find point with slope == half of what it is at max. curvature?

% -- I think this method will be way off for very high SNR signals...it needs
% -- to consider not just shape of curve, but vertical position...

% -- Possibly I shouldn't just look at individual traces, but consider 
% -- multiple events at once (at a given receiver) 







