% -- 2020-05-11
% -- Relative weight of each part of matrix (CT,CC,sP,A0)
% -- For each event... sum norms of equation coefficients?
mtrxW = zeros(Ne,4);
for ie = 1:Ne
    j1 = find(ACC(:,ie));
    j2 = find(ACT(:,ie));
    j3 = find(AsP(:,ie+2*Ne));
    j4 = find( A0(:,ie));
    mtrxW(ie,1) =  wCC*sum(vecnorm([ACC(j1,:)]'));
    mtrxW(ie,2) =  wCT*sum(vecnorm([ACT(j2,:)]'));
    mtrxW(ie,3) =  wsP*sum(vecnorm([AsP(j3,:)]'));
    mtrxW(ie,4) =   w0*sum(vecnorm([ A0(j4,:)]'));
end

% -- Ratio of triple-difference equation weight (CC and CT) over
% -- the "damping" equations
relMtrxW = sum(mtrxW(:,1:2),2)./sum(mtrxW(:,4:5),2);