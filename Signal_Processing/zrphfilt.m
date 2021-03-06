function y = zrphfilt(x,dt,lf,hf,np);
%  function y = zrphfilt(x,dt,lf,hf,np);
%
%   2020-03-27
%   Altering to take in columns, not rows!
%
%   Takes a time series sampled at DT 
%	and filters it with a 2nd order, 2 pass butterworth filter 
%   between frequencies LF and HF. If X is a matrix ZRPHFILT 
%   filters the individual rows of X. NP was added to allow different
%   numbers of poles to be used, but the default is 2.

%   It will now work as a lowpass or highpass filter if
%   lf <= 0 or hf >=nyq.

% -- Set default number of poles to 2  
if ~exist('np','var')
    np = 2;
end        

[lx,nx] = size(x);
x   = [x; zeros(lx,nx)];
nyq = 0.5/dt;

% -- Check if bandpass, lowpass, or highpass is needed
if lf > 0 && hf < nyq 
    wn    = [lf/nyq,hf/nyq];   
    [b,a] = butter(np,wn);
elseif hf < nyq  
    wn    = [hf/nyq];    
    [b,a] = butter(np,wn,'low');
elseif lf > 0   
    wn    = [lf/nyq];     
    [b,a] = butter(np,wn,'high');  
else
    y = x(1:lx,:);
    return  
end 

% -- Do filtering
y = filtfilt(b,a,x);
y = y(1:lx,:);
