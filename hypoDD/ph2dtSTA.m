function mst = ph2dtSTA(phadir,sta,maxC)
% function mst = ph2dtSTA(phadir,sta,maxC)
%
% 2020-10-21
% This function mimics the program 'ph2dt' that forms double-differences
% as input for "hypoDD"
%
%  INPUTS
%
%    phadir == [EV,STA,PHA,T,W]
%       sta == [latS,lonS]
%      maxC == maximum number of comparisons per phase
%              Note it will only cull comparisons between phases that 
%              both have more than maxC comparisons. It will keep...
%
% OUTPUTS
%
%     mst = [EV,STA1,STA2,PHA,DT,(CC)]

ue = unique(phadir(:,1));
Ne = length(ue);

M = cell(Ne,1);

for ie = 1:Ne
    
    ev = ue(ie);
    phaE = phadir(phadir(:,1)==ev,:);
    P = phaE(phaE(:,3)==1,:);
    S = phaE(phaE(:,3)==2,:);
    NP = size(P,1);
    NS = size(S,1);
    
    m = ones(NP*(NP-1)/2 + NS*(NS-1)/2, 6);
    m(:,1) = ev;
    
    jj = 0;
    for iss = 1:NP-1    
        Nj = NP-iss;
        m(jj+(1:Nj),[2 3 5]) = [repmat(P(iss,2),Nj,1) P(iss+1:end,2) ...
                                repmat(P(iss,4),Nj,1)-P(iss+1:end,4)];
        jj = jj+Nj;
    end
    
    m(jj+1:end,4) = 2;
    for iss = 1:NS-1    
        Nj = NS-iss;
        m(jj+(1:Nj),[2 3 5]) = [repmat(S(iss,2),Nj,1) S(iss+1:end,2) ...
                                repmat(S(iss,4),Nj,1)-S(iss+1:end,4)];
        jj = jj+Nj;
    end
    
    
    M{ie} = m;
end

mst = vertcat(M{:});