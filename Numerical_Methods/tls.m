function [B,sv] = tls(X,Y)
% function [B,sv] = tls(X,Y)
%
% This is a Total Least-Squares package MGB gave me...not sure where he got it. 
% 
% It's solving X*B = Y. 
%
% Can be used to get relative amplitudes between two similar vectors X,Y, as:
%  tls(X,Y). This would give the ratio A(Y)/A(X). 
%
%
[m n]   = size(X);            % n is the width of X (X is m by n)
Z       = [X Y];              % Z is X augmented with Y.
[U S V] = svd(Z,0);           % find the SVD of Z.
VXY     = V(1:n,1+n:end);     % Take the block of V consisting of the first n rows and the n+1 to last column
VYY     = V(1+n:end,1+n:end); % Take the bottom-right block of V.
B       = -VXY/VYY;
sv      = diag(S);
end
