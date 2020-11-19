function err = bootstrapHypoTD(x,Niter,AA,bb,wq,wA,minZ,maxZ)
% function err = bootstrapHypoTD(x,Niter,AA,bb,wq,wA,,minZ,maxZ) 
%
% 2020-09-16
% 
% Waldhauser talks about bootstrapping residual...
% I need to figure this out...
%
% Do I need to ignore sP???
%
% W*A*x = W*b
%
%  INPUTS
%
%     x = solution vector
% Niter = number of bootstrap iterations...
%    AA = cell array of matrices {ACC,ACT,AsP,A0,A00,AMZ}
%    bb = cell array RHS vectors {bCC,bCT,bsP,b0,b00,bMZ}
%    wq = cell array of equation weights for {ACC,ACT}
%    wA = vector of scalars [wCC,wCT,wsP,w0,w00,wMZ]
%  minZ = 
%  maxZ = 
%
%  OUTPUTS
%   
%      err = standard deviation of elements of x...
%            only meaningful for well resolved events...
%


ACC = AA{1};
ACT = AA{2};
AsP = AA{3};
A0  = AA{4};
A00 = AA{5};
AMZ = AA{6};

bCC = bb{1};
bCT = bb{2};
bsP = bb{3};
b0  = bb{4};
b00 = bb{5};
bMZ = bb{6};

wCC = wA{1};
wCT = wA{2};
wsP = wA{3};
%w0  = wA{4};
%w00 = wA{5};
w0 = 0;
w00 = 0;
wMZ = wA{6};

ACC = diag(wq{1}.^-1)*ACC;
ACT = diag(wq{2}.^-1)*ACT;
bCC = bCC.*wq{1}.^-1;
bCT = bCT.*wq{2}.^-1;

ACCx = ACC*x;
ACTx = ACT*x;
AsPx = AsP*x;

rCC = ACCx-bCC;
rCT = ACTx-bCT;
rsP = AsPx-bsP;

A  = [wCC*ACC; wCT*ACT; wsP*AsP; w0*A0; w00*A00; wMZ*AMZ];
b2 = [w0*b0; w00*b00; wMZ*bMZ];


M = length(x);
X = zeros(M,Niter);

NCC = length(bCC);
NCT = length(bCT);
NsP = length(bsP);
tic
for ii = 1:Niter

    b = [wCC*(ACCx+rCC(randi(NCC,[NCC,1])))
         wCT*(ACTx+rCT(randi(NCT,[NCT,1])))
         wsP*(AsPx+rsP(randi(NsP,[NsP,1])))
         b2];
    X(:,ii) = solveHypoTD(A,b,x,10,20,minZ,maxZ);

end
toc
err = std(X')';

% -- How to decide which events have well-resolved error estimates??


