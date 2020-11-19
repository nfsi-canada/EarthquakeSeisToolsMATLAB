% 2020-05-14
% This pice of hypoTD checks for outliers and removes them after each 
% iteration. This requires thresholds (multipliers of standard deviations)
% NsgCC,NsgCT, and NsgSP to de defined. It first reweights the equations,
% downweighting high residuals
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- CC
eCC  = WCC*(ACC*x-bCC);
sgCC = std(eCC);
jo   = find(abs(eCC) > NsgCC*sgCC);

mttcc(jcc,7) = mttcc(jcc,7).*(abs(eCC)/sgCC).^-(1/3);
mttcc(jcc(jo),:) = [];
NcullCC = NcullCC + length(jo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- CT
eCT  = WCT*(ACT*x-bCT);
sgCT = std(eCT);
jo   = find(abs(eCT) > NsgCT*sgCT);

mttct(jct,7) = mttct(jct,7).*(abs(eCT)/sgCT).^-(1/3);
mttct(jct(jo),:) = [];
NcullCT = NcullCT + length(jo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - sP
eSP  = AsP*x-bsP;
sgSP = std(eSP);
jo   = find(abs(eSP) > NsgSP*sgSP);
sP(jo,:) = [];

NcullSP = NcullSP + length(jo);