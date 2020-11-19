function [C, rep] = upgma(A, ccmin)
% function [C, rep] = upgma(A ccmin)
%
% 2020-07-17
% Performs the clustering analysis 
%    "Unweighted Pair Group Method w/ Arithmetic Mean"
% This is an agglomerative (bottom up) hierarchical clustering method
% This implentation assumes high values (near 1) are similar,
% low values (near 0) different
%
% This code was modelled after Tim Hayward's version, which he sent to me  
% on 2018-02-01
%
% Note that this assumes all the reported CCs are positive. 
%
%   INPUTS:
%
%         A == similarity matrix with similarity for each pair of events
%              (matrix of floats). Use lower triangle of a square matrix.
%              typically for me the similarity will be CC coefficients.
%              e.g. [0    0    0
%                    CC12 0    0
%                    CC13 CC23 0];
%         
%    ccmin  == threshold for including item in cluster (0 min, 1 max)
%              (single float)
%
%   OUTPUTS:
%
%        C  == cell matrix containing cluster ids.
%              Each column is one cluster.
%       rep == 'report', cell matrix containing order of clustering and 
%              correlation of each join.


% -- Get number of events
ne = size(A,1);

% -- Ensure main diagonal, upper triangle are zero
% -- Make copy of the original A
A  = tril(A,-1);
A0 = A;

% -- Initialize clusters as each event in a seperate cluster
% -- Store each cluster within it as a column vector
C = mat2cell([1:ne]',ones(ne,1),1);

% -- rep is structured [id1,id2,cc]
% -- nc is number of comparison made
rep = cell(ne*(ne-1)/2,3);
nc  = 0;

% -- Loop through CC matrix (A) and create clusters
while 1

   % -- Find initial maximum value
   [mxval, mxind] = max(A(:));
   [mxrow, mxcol] = ind2sub(size(A), mxind(1));

   % -- Stop if the cross-correlations are below threshold
   if mxval < ccmin 
        break
   end

   % -- Write comparison in report
   nc        = nc+1;
   rep(nc,:) = {C{mxrow},C{mxcol},mxval};

   % -- Get number of events in each previous cluster
   ner = length(C{mxrow});   
   nec = length(C{mxcol});

   % -- Update A with new average CCs
   j1 = 1:(mxcol-1);
   j2 = mxcol+1:mxrow-1;
   j3 = mxrow+1:size(A,1);
   A(mxcol,j1) = (nec*A(mxcol,j1)+ner*A(mxrow,j1) )/(ner+nec);
   A(j2,mxcol) = (nec*A(j2,mxcol)+ner*A(mxrow,j2)')/(ner+nec);
   A(j3,mxcol) = (nec*A(j3,mxcol)+ner*A(j3,mxrow) )/(ner+nec);

   % -- It should be obvious which
   % -- Form new A by removing 
   A(mxrow,:) = [];
   A(:,mxrow) = [];
 
   % -- Combine clusters in list
   C{mxcol} = sort([C{mxcol};C{mxrow}]);
   C(mxrow) = [];
   
   % -- Stop if we've put all events into one clusters
   if ~length(A)
        break
   end  
end

% -- Delete unused rows of report
rep(nc+1:end,:) = [];



