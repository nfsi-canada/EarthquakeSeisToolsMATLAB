% -- 2020-05-11
% -- This piece of hypoTD computes new origin times for each event, after
% -- the relative hypocenter inversion is complete
T = T0;
for ie = 1:Ne
    jp    = find(phadir(:,1)==ie);
    jss   = phadir(jp,2);
    pha   = phadir(jp,3);
    tObs  = phadir(jp,4);
    tPred = tppES(Ne*(jss-1)+ie).*vpvs.^(pha-1);
    dt0   = mean(tObs-tPred);
    T(ie) = T0(ie)+seconds(dt0);
end