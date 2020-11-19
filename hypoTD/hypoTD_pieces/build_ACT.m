% 2020-05-05
% This piece builds the matrix to invert "catalog-time" 
% triple-differences.

evA = MCT(:,1);
evB = MCT(:,2);
st1 = MCT(:,3);
st2 = MCT(:,4);

% -- Need to interpolate local velocity for each event
vA = interp1(zrng,vz,hyp(evA,3))./(vpvs.^(MCT(:,5)-1));
vB = interp1(zrng,vz,hyp(evB,3))./(vpvs.^(MCT(:,5)-1));
vA = repmat(vA,3,1);  
vB = repmat(vB,3,1);

% -- Weight equations by CC and inter-event distance
%wghtCT  = MCT(:,7)./max(dCT.^(2/3),1);
%wghtCT3 = repmat( repmat(wghtCT,3,1) ,2,1);
WCT = diag( MCT(:,7)./max(dCT.^(2/3),1) );

% -- For the equations to be valid (i.e., relate to x1-x2)
% -- me1a and me1b should be the same...
jA1 = sub2ind([Ne,Ns],evA,st1);
jA2 = sub2ind([Ne,Ns],evA,st2);
jB1 = sub2ind([Ne,Ns],evB,st1);
jB2 = sub2ind([Ne,Ns],evB,st2);

% -- Order matrix elements [gxA; gyA; gzA], and [gxB; gyB; gzB]
% -- Order matrix elements [gxA; gyA; gzA], and [gxB; gyB; gzB]
% -- Average the predicted direction cosines from the two events  
me =  0.5*( [cH(jA1).*sV(jA1); sH(jA1).*sV(jA1); cV(jA1)]        ...
           -[cH(jA2).*sV(jA2); sH(jA2).*sV(jA2); cV(jA2)] )./vA  ...
     +0.5*( [cH(jB1).*sV(jB1); sH(jB1).*sV(jB1); cV(jB1)]        ...
           -[cH(jB2).*sV(jB2); sH(jB2).*sV(jB2); cV(jB2)] )./vB;
       
j1  = repmat([1:Nct]',6,1);
j2  = [evA; evA+Ne; evA+2*Ne; evB; evB+Ne; evB+2*Ne];

ACT = sparse(j1,j2,[-me; me],Nct,3*Ne);
bCT = MCT(:,6);
%ACT = sparse(j1,j2,wghtCT3.*[-me; me],Nct,3*Ne);
%bCT = wghtCT.*MCT(:,6);

if length(find(isnan(me)))
    keyboard
end

% -- Weight equations by CC and inter-event distance
%wghtCT = MCT(:,7)./max(sqrt(dCT),1);
%wghtCT = repmat(repmat(wghtCT,3,1),2,1);
%meA = ([cH(jA1).*sV(jA1); sH(jA1).*sV(jA1); cV(jA1)] ...
%      -[cH(jA2).*sV(jA2); sH(jA2).*sV(jA2); cV(jA2)])./vA;
%meB = ([cH(jB1).*sV(jB1); sH(jB1).*sV(jB1); cV(jB1)] ...
%      -[cH(jB2).*sV(jB2); sH(jB2).*sV(jB2); cV(jB2)])./vB; 