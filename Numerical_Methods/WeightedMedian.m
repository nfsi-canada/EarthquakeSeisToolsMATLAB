function [med,xb] = WeightedMedian(x,w,wb)
% function [med,xb] = WeightedMedian(x,x,wb)
%
% 2020-07-09
% This function returns a weighted median given a data vector and its 
% corresponding weights/probabilities. It can also return error bounds
% if "percentiles" are given in wb. 
%
%    INPUTS
%
%      x = data vector
%      w = weight/probability of each element in x
%          it's critical that all weights are > 0 ...
%     wb = "bounds" of weight/probability, e.g. put [0.025,0.975] to
%          return bounds of 95% confidence interval
%
%   OUTPUTS
%
%      med = weighted median of x 
%            s.t. ha1f of weight is below, half is above
%       xb = x-position of "bounds" specified in wb. It is possible for
%            these to return NaN if e.g. cw(1) > xb(1).
%
%
%
%
%


[Nx,Nc] = size(x);
[x,jj]  = sort(x);

w = w(jj+repmat(Nx*[0:Nc-1],Nx,1));
w = w./repmat(sum(w),Nx,1);

cwL = cumsum(w);
cwU = flipud(cumsum(flipud(w)));

medL = zeros(1,Nc);
medU = zeros(1,Nc);
med  = zeros(1,Nc);

if nargin > 2
    Nb = length(wb);
    xbL = zeros(Nb,Nc);
    xbU = zeros(Nb,Nc);
else
    Nb = 0;
end

for ii = 1:Nc
    
    jL = find(cwL(:,ii)>=0.5,1);
    jU = find(cwU(:,ii)>=0.5,1,'last');

    jL1 = max(jL-1,1);
    jU1 = min(jU+1,Nx);
    
    medL(ii)  = x(jL,ii) - (cwL(jL,ii)-0.5) * diff(x([jL1,jL],ii))/diff(cwL([jL1,jL],ii));
    medU(ii)  = x(jU,ii) - (cwU(jU,ii)-0.5) * diff(x([jU,jU1],ii))/diff(cwU([jU,jU1],ii));
    
    for ib = 1:Nb
        jL = find(cwL(:,ii)>=wb(ib),1);
        jU = find(cwU(:,ii)>=1-wb(ib),1,'last');

        jL1 = max(jL-1,1);
        jU1 = min(jU+1,Nx);

        xbL(ib,ii)  = x(jL,ii) - (cwL(jL,ii)-wb(ib)) * diff(x([jL1,jL],ii))/diff(cwL([jL1,jL],ii));
        xbU(ib,ii)  = x(jU,ii) - (cwU(jU,ii)-(1-wb(ib))) * diff(x([jU,jU1],ii))/diff(cwU([jU,jU1],ii));
    end
end


medL(isnan(medL)) = x(1,find(isnan(medL)));
medU(isnan(medU)) = x(end,find(isnan(medU)));
med = (medL+medU)/2;

if Nb
    for ib = 1:Nb
        jnanL = find(isnan(xbL(ib,:)));
        jnanU = find(isnan(xbU(ib,:)));
        xbL(ib,jnanL) = x(1,  jnanL);
        xbU(ib,jnanU) = x(end,jnanU);
    end
    xb = (xbL+xbU)/2;
end

%{
% -- Ensure x is in ascending order (need not be even grid)
[X,jj] = sort(X);

Nc  = size(X,2);
med = zeros(Nc,1);
if nargin > 2
    xb = zeros(length(wb),Nc);
end

for ii = 1:Nc
    x = X(:,ii);
    w = W(jj(:,ii),ii);
    w = w/sum(w);

    % -- Eliminate points that are so low they are likely cause numerical  
    % -- errors
    jd = find(w < 2*eps);


    if length(jd) < length(x)-1

        % What if I solve for Weighted Median as a minimum of the L1 norm?
        % Start with a range around weighted mean
        %mn = sum(w.*x);

        
        x(jd) = [];
        w(jd) = [];

        % -- Need to average upper and lower weighted medians
        cwl = cumsum(w);
        cwu = flip(cumsum(flip(w)));

        % -- Use linear interpolation to find median 
        % -- (average upper and lower medians)?
        lm = interp1(cwl,x,0.5);
        um = interp1(cwu,x,0.5);
        lm(isnan(lm)) = min(x);
        um(isnan(um)) = max(x);

        med(ii) = (lm+um)/2;

        % -- Return other percentiles of weights
        if nargin > 2
            lb = interp1(cwl,x,wb);
            ub = interp1(cwu,x,1-wb);
            lb(isnan(lb)) = min(x);
            ub(isnan(ub)) = max(x);
            xb(:,ii) = (lb+ub)/2;
        end

    else

        % -- For extreme cases there is only one effective non-zero point
        % -- Assume curve intersects zero-halfway between median and 
        % -- neighboring points
        [~,jm] = max(w); 
        med(ii) = x(jm);
        if nargin > 2
            xf = ([x(jm-1) x(jm+1)]+med)/2;
            xb(:,ii) = interp1([0 1],xf,wb,'pchip');
        end
    end      
end
%}
