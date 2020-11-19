function y = pad0(x,n)
% function y = pad0(x,n)
% 
% 2020-07-17
% This function just zero-pads columns of an input matrix "x" such that
% it has "n" rows. If "x" is already taller than "n" rows, it will be
% truncated to have only "n" rows.

[nr,nc]=size(x);
if(nr > n)
    y = x(1:n,:);
else
    y = [x; zeros(n-nr,nc)];
end

