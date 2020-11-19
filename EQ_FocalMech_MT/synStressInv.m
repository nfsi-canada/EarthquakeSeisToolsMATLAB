% 2020-08-06
% This script creates a set of focal mechanisms consistent with a given
% stress tensor so that the different stress inversion codes can be tested.

clear all


% -- Create stress tensor 'X' or in 5x1 vector form 'x'
% -- Need x0 = [x11 x12 x13 x22 x23]
% -- Start with sigma1,2,3 along x/y/z axes. They must sum to zero.
x123 = [-1 0 1];
X0 = diag(x123-mean(x123));

% -- Now rotate the stress tensor (start with a unit quaternion)
q  = [0 0 0 1];
q  = q/norm(q);
RM = Quat2RotMtrx(q);
X0 = RM'*X0*RM;
x0 = X0([1 2 3 5 6]');

% -- Compute actual trend/plunge (T0,P0); and shape ratio R0
% -- Nnote that sigma1 (most compressive) is most negative
[V,S] = eig(X0);
[S,jS] = sort(diag(S),'descend');
V      = V(:,jS);
V(:,V(3,:)>0) = -V(:,V(3,:)>0);
P0 = -asind(V(3,:))';

T0 = 90-atan2d(V(2,:),V(1,:))';
T0(T0<0) = T0(T0<0)+360;
R0 = (S(1)-S(2))/(S(1)-S(3)); 

% -- Select number of events (N), create random fault planes...
% -- The normals of the fault planes are what we need to build the
% -- matrix A... (normals pointing up)
N = 20;
n = zeros(N,3);
n00 = [0 0 1];
for ii = 1:N
    R = Quat2RotMtrx(RandomQuaternion(1));
    n(ii,:) = n00*R;
end

% -- Ensure normals point up
n(n(:,3)<0,:) = -n(n(:,3)<0,:);

% -- Build stress inverse matrix...
A = BuildStressMatrix(n);

% -- Compute rakes (in RHS vector b)
v = reshape(-A*x0,[N 3]);
v = v./repmat(vecnorm(v')',1,3);

% -- Convert normal,slip vectors into strike/dip/rake
sdr = nv2sdr(n,v);

% -- THIS CAN'T BE SOLVED PERFECTLY!!
% -- Angelier [2002,GJI] has some complicated method that might account
% -- for all of this...I don't understand it all
[TRND,PLNG,R,S,msft] = StressInvM1984(sdr)

%[TRND,PLNG,R,uS,uR] = StressInvV2014(sdr,100,0.1,'C',[68 95])


%{
x = AngelierStressSolver(n,v,randn(5,1),0.11);

X = x([1 2 3; 2 4 5; 3 5 1]);
X(3,3) = -X(1,1)-X(2,2);
[V,S]  = eig(X);
[S,jS] = sort(diag(S));
V = V(:,jS);
V(:,V(3,:)>0) = -V(:,V(3,:)>0);   
R = (S(1)-S(2))/(S(1)-S(3));  

% -- Compute plunge and trend, ensure 0<trend<360 
PLNG = -asind(V(3,:))';
TRND = 90-atan2d(V(2,:),V(1,:))';
TRND(TRND<0) = TRND(TRND<0)+360;
[TRND,PLNG]
R
%}