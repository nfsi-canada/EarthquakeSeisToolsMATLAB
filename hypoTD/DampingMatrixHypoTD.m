function [GD,dD] = DampingMatrixHypoTD(C,Ne0,lmbd);
% function [GD,dD] = DampingMatrixHypoTD(C,Ne0,lmbd);
% 
% 2021-01-15
% This function produces a damping matrix for double-difference hypocenter
% relocation, with weights tailored to the size of individual clusters.
%
% This was modified from "DampingMatrixHypoTD.m" to remove time columns
%
% INPUTS
%
%     C  == [Ne x 1] vector of cluster IDs 
%    Ne0 == original number of events
%     W0 == Weight of zero-sum constraits
%
% OUTPUTS
%
%     GD == damping matrix
%     dD == all-zero vector
%
Ne = length(C);
uC = unique(C);
NC = length(uC);

gd = sparse(Ne,1);
for ii = 1:NC
    jj = find(C == uC(ii));
    gd(jj) = lmbd*length(jj)/Ne0;
end
GD = diag(repmat(gd,3,1));
dD = sparse(3*Ne,1);