function sfn = shift(fn,dt,t0)
% function sfn = shift(fn,dt,t0)
%
% 2020-03-27
% This is updated from MGB's function to apply to columns instead of 
% rows.
%
% FUNCTION [SFN] = SHIFT(FN,DT,T0)
% Function SHIFT produces (circularly) shifted time series SFN from
% input time series FN where T0 is the desired time
% shift, DT is the sampling interval. If array of 
% traces is passed then T0 should be an array too.
% Uses FFT (ie sinc(x)) interpolation.
%
%

% -- Preliminaries.
[N,M] = size(fn);
w = (2*pi/(N*dt))*[0:N-1]';

% -- If only one time shift was provided (impyling a constant shift for
% -- each column), replicate it my the M columns
if length(t0) == 1
    t0 = t0*ones(M,1);
end

% -- Take fourier transform
ffn = fft(fn,N);

% -- Loop through traces. Check for even/odd length since it affects
% -- definition of Nyquist frequency.
if mod(N,2) == 0
    for ii = 1:M
        ffn(1:N/2,ii)   = ffn(1:N/2,ii).*exp(-i*w(1:N/2)*t0(ii));
        ffn(N/2+1,ii)   = ffn(N/2+1,ii)*cos(pi*t0(ii)/dt);
        ffn(N/2+2:N,ii) = conj(ffn(N/2:-1:2,ii));
    end    
else
    for ii=1:M
        ffn(1:(N+1)/2,ii)   = ffn(1:(N+1)/2,ii).*exp(-i*w(1:(N+1)/2)*t0(ii));
        ffn((N+1)/2+1:N,ii) = conj(ffn((N+1)/2:-1:2,ii));
    end
end

dum = real(ifft(ffn,N));
sfn = dum(1:N,:);