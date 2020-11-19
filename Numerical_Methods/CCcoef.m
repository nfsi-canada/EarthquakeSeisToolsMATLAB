function cc = CCcoef(x,y)
% function cc = CCcoef(x,y)
%
% 2016-09-22
% This function finds a simple correlation coefficient between
% two equal-length vectors. If matrices are provided,
% a CC is produced for each column (returned in a row vector).

cc = sum(x.*y) ./ sqrt(sum(x.^2).*sum(y.^2));

