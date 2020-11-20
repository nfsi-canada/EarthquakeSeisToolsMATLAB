function [sdr1,sdr2,err,NR] = FocMechAPP(PFM,PVH,crsX,xyzP,xyzI,xyzT,Qp,Qs,vpvs)
% function [sdr1,sdr2,err,NR] = FocMechAPP(PFM,PVH,crsX,xyzP,xyzI,xyzT,Qp,Qs,vpvs)
%
% 2020-11-19
% This function estimates an earthquake focal mechanism from
% P-wave first-motions and absolute-value P-SV-SH amplitude ratios. It 
% requires a pre-defined grid that uniformly samples all possible focal
% mechanisms, as it performs a grid search. 
%
% See "createDCgrid.m" to build a uniform grid of double-couple (DC)
% moment tensors. 
%
%   INPUTS
%
%    PFM = [AZI,TOA,+1/-1,SNR] (see column descriptions below)
%    PVH = [AZI,TOA,TP,TS,FRQ,Ap,Av,Ah] (see column descriptions below)
%   crsX = 6 x Ng list of grid-search (DC) MTs [M11 M22 M33 M12 M13 M23]'
%   xyzP = 3 x Ng matrix of pressure axes corresponding to [x y z]
%   xyzI = 3 x Ng matrix of intermediate axes corresponding to [x y z]
%   xyzT = 3 x Ng matrix of tension axes corresponding to [x y z]
%     Qp = P-wave attenuation (unitless), default = 1350
%     Qs = S-wave attenuation (unitless), default = 1350
%   vpvs = vp/vs ration, default 1.75
%
%    column descriptions:
%
%        AZI == azimuth, degrees east of north
%        TOA == take-off angle in degrees
%               0 == up 
%               90 == horizontal
%               180 == down
%     +1/-1 == either +1 for up-first motion or -1 for down-first
%       SNR == signal-to-noise ratio
%        TP == P travel time (s)
%        TS == S travel time (s)
%       FRQ == Nominal frequency of bandpass amplitudes were 
%              measured at (I recommend using the highpass frq.)
%        Ap == P absolute amplitude  
%        Av == SV absolute amplitude
%        Ah == SH absolute amplitude (or zero for single-component station)
%
%  OUTPUTS
%
%    sdr1,sdr2 = [strike,dip,rake] of best-fit focal mechanism
%          err = standard error, in degrees
%                i.e. the true mechanism is on average "err" degrees away
%                from the reported one with this degree of uncertainty
%           NR = total number of amplitude ratios



% -- Check if optional parameters were provided
if nargin < 9
    vpvs = 1.75;
    if nargin < 8
        Qs = 600;
        if nargin < 7
            Qp = 1350;
        end
    end
end

% -- nx = 6 == Number of independent MT elements [M11 M22 M33 M12 M13 M23]
[nx,ng] = size(crsX);

% -- Build first motion linear system
[AFM,bFM,NFM] = build_AFM(PFM);

% -- Convert first-motions fits to probability
% -- Note the somewhat arbitrary parameters turning SNR into a weight
wFM    = 0.9*tanh(PFM(:,4)-0.5);
swFM   = 0.5*abs(sign(AFM*crsX)-repmat(bFM,1,ng));
ProbFM = 1 - swFM.*repmat(wFM,1,ng).*abs(AFM*crsX).^(1/2);
ProbFM = prod(ProbFM);

% -- Build P-SV-SH absolute amplitude ratio linear system
[AR,bR,NR] = build_AR(PVH,Qp,Qs,vpvs);

% -- Convert linear system misfits to (relative) probability
msftR = vecnorm( abs((AR(:,1:nx)*crsX)./(AR(:,nx+1:end)*crsX))-repmat(bR,1,ng),1);
msftR = msftR/min(msftR);     
ProbR = msftR.^(-(NR));


% -- 
ProbX  = ProbFM.*ProbR;
ProbX  = ProbX/sum(ProbX);
        
% -- Find maximum probability grid-point
% -- Take Pressure, Tension, 
[mxPrb,jmax] = max(ProbX);
X0 = [xyzP(:,jmax) xyzT(:,jmax) xyzI(:,jmax)];

% -- For efficiency ignore points with very low probabilities
jj    = find(ProbX > mxPrb*1e-6);
ProbX = ProbX/sum(ProbX(jj));

% -- Search around best grid 
[sdr1,sdr2,err] = solveDC_L1(X0,ProbX(jj),xyzP(:,jj),xyzT(:,jj),xyzI(:,jj));


