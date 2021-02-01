function [GZ,dZ] = ZeroShiftHypoTD(C,W0,evsSP)
% function [GZ,dZ] = ZeroShiftHypoTD(C,W0,evsSP)
%
% 2020-01-15
% This function forms a matrix and vector to enforce that the mean X/Y/Z
% values all remain constant, for double-difference hypocenter inversion.
% It takes in a list of ClusterIDs, to enforce the zero mean constraints
% separately on each cluster.
%
% This was modified from "ZeroShiftHypoDD.m"
%
% INPUTS
%
%     C  == [Ne x 1] vector of cluster IDs 
%     W0 == Weight of zero-sum constraits
%  evsSP == (optional) event indices with sP constraints
%    
% OUTPUTS
%
%     GZ == matrix with a zero-sum constrait for each dimension, for
%           each cluster
%     dZ == zeros




Ne = length(C);
uC = unique(C);
NC = length(uC);

j1 = cell(NC,1);
j2 = cell(NC,1);

% -- If no sP phases are given then all zero-sum constraints are highly
% -- weighted
if nargin < 3
    for ii = 1:NC
        jj = find(C == uC(ii));
        Nj = length(jj);
        j1{ii} = 3*(ii-1)+[ones(Nj,1); 2*ones(Nj,1); 3*ones(Nj,1)];
        j2{ii} = [jj; Ne+jj; 2*Ne+jj];
    end
    me = W0*ones(3*Ne,1);
    
% -- If clusters include an sP depth constraint put only a tiny weight
% -- on the vertical zero-sum constraint
else
    Csp  = unique(C(evsSP));
    Cn   = setdiff(uC,Csp);
    NCsp = length(Csp);
    NCn  = length(Cn);
    me = cell(NC,1);
    for ii = 1:NCsp
        jj = find(C == Csp(ii));
        Nj = length(jj);
        j1{ii} = 3*(ii-1)+[ones(Nj,1); 2*ones(Nj,1); 3*ones(Nj,1)];
        j2{ii} = [jj; Ne+jj; 2*Ne+jj];
        me{ii} = min(0.01,W0/1000)*ones(3*Nj,1);
    end
    for ii = 1:NCn
        jj = find(C == Cn(ii));
        Nj = length(jj);
        j1{ii} = 3*(ii-1)+[ones(Nj,1); 2*ones(Nj,1); 3*ones(Nj,1)];
        j2{ii} = [jj; Ne+jj; 2*Ne+jj];
        me{ii} = W0*ones(3*Nj,1);
    end
    me = vertcat(me{:});
end
    
GZ = sparse(vertcat(j1{:}),vertcat(j2{:}),me,3*NC,3*Ne);
dZ = sparse(3*NC,1);

