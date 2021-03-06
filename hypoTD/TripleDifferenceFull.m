function B = TripleDifferenceFull(A)
% function B = TripleDifferenceFull(A)
%
% 2020-05-05
% This function takes in a list of differential travel-times and computes
% triple-differences DDT = ((T1a-T1b)-(T2a-T2b))_obs, where 1,2 are 
% stations, a,b are events (and the times are travel times)
% This version does not cull the TDs to reduce redundancy...
%
% INPUTS
%
%     A = [EVa,EVb,STA,PHA,DT,CC]
%
% OUTPUTS
%
%     B = [EVa,EVb,STA1,STA2,PHA,DDT,CC]
%

if ~length(A)
    B = zeros(0,7);
    return
end
    
B = zeros(1e7,7);
jj  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- First deal with P-waves
M = A(A(:,4)==1,:);

upair = unique(M(:,1:2),'rows');

for evAB = upair'
    
    Mab  = M(find((M(:,1)==evAB(1)).*(M(:,2)==evAB(2))),:);
    Ns   = size(Mab,1);
    
    if Ns > 1
           
        jss  = nchoosek(1:Ns,2);
        ddt  = Mab(jss(:,1),5)-Mab(jss(:,2),5);
        CC   = min( Mab(jss(:,1),6),Mab(jss(:,2),6) );
        N    = size(jss,1);
        
        B(jj+[1:N],:) = [repmat(evAB',N,1),Mab(jss(:,1),3),Mab(jss(:,2),3),ones(N,1),ddt,CC];
        jj = jj + N;
       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Repeat for S-waves
M = A(A(:,4)==2,:);

upair = unique(M(:,1:2),'rows');

for evAB = upair'
    
    Mab  = M(find((M(:,1)==evAB(1)).*(M(:,2)==evAB(2))),:);
    Ns   = size(Mab,1);
    if Ns > 1
        
        jss  = nchoosek(1:Ns,2);
        ddt  = Mab(jss(:,1),5)-Mab(jss(:,2),5);
        CC   = min( Mab(jss(:,1),6),Mab(jss(:,2),6) );
        N    = size(jss,1);

        B(jj+[1:N],:) = [repmat(evAB',N,1),Mab(jss(:,1),3),Mab(jss(:,2),3),2*ones(N,1),ddt,CC];
        jj = jj + N;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Combine and sort

B(jj+1:end,:) = [];
B = sortrows(B);