function err = errorHypoTD(A,x,b)
% function err = errorHypoTD(A,x,b)
%
% 2020-07-29
% This function estimates error on the x,y,z coordinates of each 
% earthquake in a HypoTD problem, following the procedure outlined in
% Waldhauser and Ellsworth [2000,BSSA] equations 12 and 13. This is 
% only computationally feasible for small systems...
%
% INPUTS
%
%    A == final, hypoTD matrix (only CC,CT,and sP parts?)
%    x == final hypocenter estimates [x; y; z]
%    b == final hypoTD RHS vector
%
% OUTPUTS
%
%    err == estimated error of each coordinate (in whatever the input
%           unit is, typically km)
keyboard
tic
[M,N] = size(A);
N0 = N/3;

%if N > M
%    disp('System is underdetermined!!')
%    err = Inf*ones(M,1);
%end

% -- Count the number of equations for each event...
% -- Note that columns are sorted [all X, all Y, all Z]
% -- Use Z columns because sP equations only depend on Z
Nq = zeros(N0,1);
for ii = 1:N0
    Nq(ii) = nnz(A(:,(2*N0)+ii));
end

% -- Sort events/columns being kept from those being removed
Er = find(Nq < 5);
Ek = find(Nq >= 5);
jcr = [Er; Er+N0; Er+2*N0];
jck = [Ek; Ek+N0; Ek+2*N0];

% -- Shift columns over of event being removed over to the RHS?
%b = b-A(:,jcr)*x(jcr);
%A(:,jcr) = [];
%x(jcr)   = [];

Nnz = zeros(M,1);
for ii = 1:M
    Nnz(ii) = nnz(A(ii,:));
end
A = A(find(Nnz),:);
b = b(find(Nnz));
[M,N] = size(A);
N = N/3;


[U,S,V] = svd(full(A));

Sn2 = diag(diag(S).^-2);
C = diag(V*Sn2*V');

%Sn2 = S;
%Sn2(find(Sn2)) = Sn2(find(Sn2)).^-2;

% -- Build equation 13
d  = A*x-b;
md = mean(d);

%dk = d(jck);
%mdk = mean(dk);

% -- Possibly need to remove weights????
var = ( sum((d-md).^2) - sum(d-md)^2/M )/(M-3*N);
%var = ( sum((dk-mdk).^2) - sum(dk-mdk)^2/Mk )/(Mk-3*Nk);

% -- Estimate errors (equaiton 12)
%e = Inf*ones(3*N0,1);
%e(jck) = sqrt(C*var);

e = sqrt(C(jck)*var);

toc

figure(1)
clf
hold on
plot(e(1:Nk),'bo-')
plot(e([1:Nk]+Nk),'ro-')
plot(e([1:Nk]+2*Nk),'go-')
ylim([0 500])