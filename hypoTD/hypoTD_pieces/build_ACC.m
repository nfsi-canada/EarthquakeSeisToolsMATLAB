% 2020-05-05
% This piece builds the matrix to invert cross-correlation 
% triple-differences.
%
%    ((T1a-T2a)-(T1b-T2b))_obs - ((T1a-T2a)-(T1b-T2b))_cal
% == ((T1a-T1b)-(T2a-T2b))_obs - ((T1a-T1b)-(T2a-T2b))_cal
%   (stations 1,2, events a,b)

evA = MCC(:,1);
evB = MCC(:,2);
st1 = MCC(:,3);
st2 = MCC(:,4);

% -- Need to interpolate local velocity for each event
vA = interp1(zrng,vz,hyp(evA,3))./(vpvs.^(MCC(:,5)-1));
vB = interp1(zrng,vz,hyp(evB,3))./(vpvs.^(MCC(:,5)-1));
vA = repmat(vA,3,1);  
vB = repmat(vB,3,1);

% -- Weight equations by CC and inter-event distance
%wghtCC  = MCC(:,7)./max(dCC.^(2/3),1);
%wghtCC3 = repmat( repmat(wghtCC,3,1) ,2,1);
WCC = diag( MCC(:,7)./max(dCC.^(2/3),1) );

% -- For the equations to be valid (i.e., relate to x1-x2)
% -- me1a and me1b should be the same...
jA1 = sub2ind([Ne,Ns],evA,st1);
jA2 = sub2ind([Ne,Ns],evA,st2);
jB1 = sub2ind([Ne,Ns],evB,st1);
jB2 = sub2ind([Ne,Ns],evB,st2);

% -- Order matrix elements [gxA; gyA; gzA], and [gxB; gyB; gzB]
% -- Average the predicted direction cosines from the two events  
me =  0.5*( [cH(jA1).*sV(jA1); sH(jA1).*sV(jA1); cV(jA1)]        ...
           -[cH(jA2).*sV(jA2); sH(jA2).*sV(jA2); cV(jA2)] )./vA  ...
     +0.5*( [cH(jB1).*sV(jB1); sH(jB1).*sV(jB1); cV(jB1)]        ...
           -[cH(jB2).*sV(jB2); sH(jB2).*sV(jB2); cV(jB2)] )./vB;
       
% -- Order of matrix elements is [gxA; gyA; gzA], and [gxB; gyB; gzB]
j1  = repmat([1:Ncc]',6,1);
j2  = [evA; evA+Ne; evA+2*Ne; evB; evB+Ne; evB+2*Ne];

ACC = sparse(j1,j2,[-me; me],Ncc,3*Ne);
bCC = MCC(:,6);

%ACC = sparse(j1,j2,wghtCC3.*[-me; me],Ncc,3*Ne);
%bCC = wghtCC.*MCC(:,6);


%meA = ([sH(jA1).*sV(jA1); cH(jA1).*sV(jA1); cV(jA1)] ...
%      -[sH(jA2).*sV(jA2); cH(jA2).*sV(jA2); cV(jA2)])./vA;
%meB = ([sH(jB1).*sV(jB1); cH(jB1).*sV(jB1); cV(jB1)] ...
%      -[sH(jB2).*sV(jB2); cH(jB2).*sV(jB2); cV(jB2)])./vB; 
%meB=meA;
  
%meA = ([cH(jA1).*sV(jA1); sH(jA1).*sV(jA1); cV(jA1)] ...
%      -[cH(jA2).*sV(jA2); sH(jA2).*sV(jA2); cV(jA2)])./vA;
%meB = ([cH(jB1).*sV(jB1); sH(jB1).*sV(jB1); cV(jB1)] ...
%      -[cH(jB2).*sV(jB2); sH(jB2).*sV(jB2); cV(jB2)])./vB;   
%meB = meA;
%{
me1a = -[cH(jA1).*sV(jA1);
         sH(jA1).*sV(jA1);
         cV(jA1)]./vA;
me2a =  [cH(jA2).*sV(jA2);
         sH(jA2).*sV(jA2);
         cV(jA2)]./vA;
me1b =  [cH(jB1).*sV(jB1);
         sH(jB1).*sV(jB1);
         cV(jB1)]./vB;    
me2b = -[cH(jB2).*sV(jB2);
         sH(jB2).*sV(jB2);
         cV(jB2)]./vB;     
%}


%j1  = repmat([1:Ncc],3,1);
%j1  = repmat(j1(:),2,1);
%j2  = repmat([evA',evB'],3,1);
%j2  = [evA',evB'];
%j2  = [j2; j2+Ne; j2+2*Ne];
%j2  = [evA; evA+Ne; evA+2*Ne; evB; evB+Ne; evB+2*Ne];
%j2A  = [evA, evA+Ne, evA+2*Ne]'; 
%j2B  = [evB, evB+Ne, evB+2*Ne]';
%j2   = [j2A(:); j2B(:)];
%ACC = sparse(j1,j2,wghtCC.*[me1a+me2a;me1b+me2b],Ncc,3*Ne);