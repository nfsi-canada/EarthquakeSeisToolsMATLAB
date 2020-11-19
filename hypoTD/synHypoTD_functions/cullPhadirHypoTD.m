function [phadir,jkeep] = cullPhadir(phadir,NF)
% function [phadir,jkeep] = cullPhadir(phadir,NF)
%
% 2019-11-04
% This function removes a fraction of phases from phadir, with specific 
% fraction for [P,S] of phases to be removed indicated in NF;
% e.g. NF = [0.2, 0.4] will select 20% of P phases to delete and 40% of S. 
%
% phadir == [EV,STA,PHA,T] ... an N_phase x 4 matrix where PHA = 1,2 for P,S
% jkeep  == rows of phadir that were kept;

%np = size(phadir,1);

jp = find(phadir(:,3)==1);
js = find(phadir(:,3)==2);

np = length(jp);
ns = length(js);

ndp = fix(NF(1)*np);
nds = fix(NF(2)*ns);

jdp = jp(randperm(np,ndp));
jds = js(randperm(ns,nds));

phadir([jdp; jds],:) = [];

jkeep = setdiff(1:(np+ns),[jdp; jds]);