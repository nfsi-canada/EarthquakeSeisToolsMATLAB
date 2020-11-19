function tt = Kushnir1990(x)
% function tt = Kushnir1990(x)
%
% 2020-01-07
% This function is based on Kushnir et al. [1990,BSSA] autoregressive
% picking technique for a single trace...
%
% Follow Pisarenko et al. [1987,PEPI] procedure...between eq. 9 and 10?

% -- Degree of autoregressive model
n = 2;

[N,D] = size(x);

L = zeros(N,1);
x = x/max(abs(x(:)));

% -- Loop over time samples
for ii = n+1:N-n-1
        
    x1 = x(1:ii,:);
    x2 = x(ii+1:end,:);
    
    % -- Cross-correlate, only need up to a lag of 'n'
    R1 = xcorr(x1,n);
    R2 = xcorr(x2,n);
    R1(1:n,:) = [];
    R2(1:n,:) = [];
    
    % -- Apply Levinson-Durbin procedure to the the autocovariances?
    A1 = levinson(R1,n)';
    A2 = levinson(R2,n)';

    % -- Equations just before (10) in Pisarenko et al.
    % -- reshape into 3x3 matrix in the same line (if there are 3
    % -- components!
    B1 = reshape(R1(1,:)-dot(R1(2:end,:),A1(2:end,:))/ii    ,D,D);
    B2 = reshape(R2(1,:)-dot(R2(2:end,:),A2(2:end,:))/(N-ii),D,D);
    
    % -- Compute likelihood (divide exponents by N for numerical stability)
    L(ii) = det(B1)^(-(ii)/N) * det(B2)^(-(N-ii)/N);
      
end

jj     = [n+1:N-n-1]';
L      = real(L);
L(jj)  = L(jj)-min(L(jj));
[~,tt] = max(L);


%tt = jt;

%{
q = zeros(N,1);
q(jj) = -(jj-(n+1)).*(jj-(N-n-1));
q = q*sum(L)/sum(q);
Lq = L;
Lq(jj(2:end-1)) = max(L(jj(2:end-1))./q(jj(2:end-1))-1,0);

% -- I want some measure of the 'spread'
% -- This is akin to a standard deviation
% -- Should I change it to L.^2 ??
sig = sqrt( sum(Lq.^2.*([1:N]'*dt-tt).^2 )/sum(Lq.^2))


figure(1)
clf
subplot(211)
hold on
plot([1:N]*dt,L)
plot([tt,tt],[0,max(L)],'r--')
plot([1:N]*dt,q,'g')
subplot(212)
hold on
plot([1:N]*dt,x)
plot([tt,tt],[-1,1],'r--')
keyboard
%}
% -- Construct predicted timeseries with toeplitz matrix
%B1 = zeros(D);
%B2 = zeros(D);
%for j1 = 1:D
%    for j2 = 1:D 
%        kk = 3*(j1-1)+j2;
%        
%        B1(j1,j2) = (R1(1,kk)-dot(R1(2:n+1,kk),A1(kk,2:n+1)))./(ii);
%        B2(j1,j2) = (R2(1,kk)-dot(R2(2:n+1,kk),A2(kk,2:n+1)))./(N-ii);
%    end
%end

