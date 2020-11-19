function [trnd,plng,R,S,msft] = StressInvM1984(sdr)
% function [trnd,plng,R,S,msft] = StressInvM1984(sdr)
%
% 2019-11-27
% This function does a linear inversion to estimate a stress tensor given
% a list of earthquake focal mechanisms (strike,dip,rake), following the
% method of Michael [1984,JGR].
%
% The stress axes (sigma1,sigma2,sigma3) are ordered from most compressive
% to least compressive, with a tension==positive sign convention 
%
%    INPUTS
%
%   sdr = an NE x 3 matrix of [STRIKE,DIP,RAKE] in degrees
%         strike: 0 to 360 in degrees E of N
%         dip:    0 to 90 (right-hand rule from strike)
%         rake:   -180 to +180, where 0 is dextral strike-slip
%                                     +90 is reverse
%                                     -90 is normal
%
%    OUTPUTS
%
%   trnd = 'trend' of sigma1-3, in degrees E of N
%   plng = 'plunge' of sigma1-3, in degrees E of N
%      R = shape ratio (sig1-sig2)/(sig1-sig3)
%              R=0: uniaxial tension
%              R=1: uniaxial compression
%              NOTE: R = 1-phi, (M1984 report phi)
%      S = principal stresses [sig1,sig2,sig3] (compressional--tensional)
%   msft = L2-misfit  norm(A*x-b)/norm(b)


% -- Get fault normals (n) and slip vectors (v) in ENU convetion
% -- Normals point up, slip vector describes motion of hanging wall
[n,v] = sdr2nv(sdr);

% -- Construct A with Vavrycuk's [2014,GJI,Eq 8] notation (but not transposed).
% -- Row order is all E/s1 equations, then N/s2, then U/s3
% -- Column order is [S11, S12, S13, S22, S23]
A = BuildStressMatrix(n);
b = -v(:);

% -- Solve for stress tensor, reform into 3x3   
disp('USING AN ITERATIVE STRESS SOLVER!')
disp('Not simply a straightforward linear inversion')
[x,msft] = iterativeStressSolver(A,b,1e-6);

%x = A\b;
X = x([1 2 3; 2 4 5; 3 5 1]);
X(3,3) = -X(1,1)-X(2,2);

% -- Decompose the stress tensor, ensure eigenvalues go from most 
% -- negative/compressional to most positive/tensional. 
% -- Ensure eigenvetors point down (i.e. that their 3rd element is negative). 
[V,S]  = eig(X);
S = flip(diag(S));
V = fliplr(V);
V(:,V(3,:)>0) = -V(:,V(3,:)>0);   
R = (S(1)-S(2))/(S(1)-S(3));  

% -- Compute plunge and trend, ensure 0<trend<360 
plng = -asind(V(3,:))';
trnd = 90-atan2d(V(2,:),V(1,:))';
trnd(trnd<0) = trnd(trnd<0)+360;

% -- L2 misfit (normalize predicted rakes?)
%msft = norm(A*x-b)/norm(b);    

%bp = reshape(A*x,[size(sdr,1),3]);
%bp = bp./repmat(vecnorm(bp')',1,3);
%bp = bp(:);
%msft = norm(bp-b)/norm(b);

% -- Closer to M1984 version
%A = [n1.*(1-n1.^2+n3.^2), n2.*(1-2*n1.^2), n3.*(1-2*n1.^2), n1.*(n3.^2-n2.^2),  -2*n1.*n2.*n3,
%     n2.*(n3.^2-n1.^2),   n1.*(1-2*n2.^2), -2*n1.*n2.*n3,   n2.*(1-n2.^2+n3.^2), n3.*(1-2*n2.^2),
%     n3.*(n3.^2-n1.^2-1), -2*n1.*n2.*n3,   n1.*(1-2*n3.^2), n3.*(n3.^2-n2.^2-1), n2.*(1-2*n3.^2)];

% -- Exactly as Vavrycuk [2014,GJI] writes it
%A = [n1.*(n2.^2+2*n3.^2),  n2.*(-n1.^2+n3.^2), n3.*(-2*n1.^2-n2.^2); ...
%         n2.*(1-2*n1.^2),     n1.*(1-2*n2.^2),        -2*n1.*n2.*n3; ...
%         n3.*(1-2*n1.^2),       -2*n1.*n2.*n3,      n1.*(1-2*n3.^2); ...
%      n1.*(-n2.^2+n3.^2), n2.*(n1.^2+2*n3.^2), n3.*(-n1.^2-2*n2.^2); ...
%           -2*n1.*n2.*n3,     n3.*(1-2*n2.^2),      n2.*(1-2*n3.^2)]'; 




