function y = bessel_filt(x,dt,lf,hf,np)
% function y = bessel_filt(x,dt,lf,hf,np)
%
% 2020-01-11
% This function perfoms Bessel fitlering as recommended by Steve Roecker
% ...
%

% -- Set default number of poles to 2  
if nargin < 5
    np = 2;
end        

L = size(x,1);
x = [x; zeros(size(x))];
F = 1/dt;
nyq = F/2;

% -- Check if bandpass, lowpass, or highpass is needed
if lf > 0 && hf < nyq 
    [z,p,k]    = besself(np,[lf,hf]*2*pi);       
    [zd,pd,kd] = bilinear(z,p,k,F);    
    sos        = zp2sos(zd,pd,kd); 
elseif hf < nyq  
    [z,p,k]    = besself(np,hf*2*pi,'low');       
    [zd,pd,kd] = bilinear(z,p,k,F);    
    sos        = zp2sos(zd,pd,kd); 
elseif lf > 0   
    [z,p,k]    = besself(np,lf*2*pi,'high');       
    [zd,pd,kd] = bilinear(z,p,k,F);    
    sos        = zp2sos(zd,pd,kd); 
else
    y = x(1:L,:);
    return  
end  

y = real(sosfilt(sos,x));
y = y(1:L,:);
