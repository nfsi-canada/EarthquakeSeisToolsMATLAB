function [DT,CC,pks] = L1align10(xA,xB,dt,win,mxlg,nmax,phase,Ndt)
% function [DT,CC,pks] = L1align10(xA,xB,dt,win,mxlg,nmax,phase,Ndt)
%
% 2019-01-31
% This code was modified from 'superCC10.m'. It computes a L1 misfit between
% waveforms that have been normalized, and corrected for polarity. Although it
% finds minima of the L1 misfit, it actually just records CCs.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018-04-25
% This function aligns waveforms that have multiple channels, finding the 
% highest absolute value CC. The name refers to the 'supertrace' composed of the 
% multiple components. In addition to finding the very highest, it also
%
% 2018-08-09
% Added a phase option. Default is for P-wave where all components should share
% a polarity. S-wave option checks the better polarity for each component---
% output CCs will always be positive. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    INPUTS
%
%       xA,xB == the data, each channel in a separate column
%                for standard use, these will be aligned such that the S-pick
%                lies on the first sample of 'win'
%          dt == sample interval in seconds 
%         win == indices of the window to use
%                must have win(1) > mxlag, win(2) < size(x,1)-mxlag
%        mxlg == maximum lag (+\-) in samples
%        nmax == number of maxima to report in 
%       phase == 'P' (default) for simple supertrace CC, but looks for high
%                |CC|. 'S' uses positive polarity for each channel. Any other 
%                input does a simple CC and looks for high positive values. 
%         Ndt == Test every 'Ndt' samples (default is 1).  
%
%    OUTPUTS
%
%         DT  == optimal tA-tB in seconds
%         CC  == correlation coefficient at optimal delay    
%        pks  == up to nmax peaks cross-correlation function [T,CC]
%
% NOTE (2018-11-19) I might be able to do this without the loop to build MA,MB.

if nargin < 8
    Ndt = 1;
    if nargin < 7
        phase = 'P';
    end
end

lg  = [-mxlg:Ndt:mxlg-Ndt];
nlg = length(lg);
%nlg = 2*mxlag;


[len,ncha] = size(xA);
lwin = length(win);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Compute CCF

% -- Build matrices so that whole CCF can be computed in one line of code.
MA = zeros(ncha*lwin,nlg);
MB = zeros(ncha*lwin,nlg);

% -- What if I didn't need the above loop?
jA = fix(lg/2);
jB = jA-lg;

for ii = 1:nlg

    % -- 'Centered' times for given lags
    %jA = fix(lg(ii)/2);
    %jB = jA-lg(ii);

    % -- v for 'vector', put different channels into one vector
    vA = xA(win+jA(ii),:);
    vB = xB(win+jB(ii),:);
    MA(:,ii) = vA(:);
    MB(:,ii) = vB(:);

end



% -- Normalize the power in each column to 1
MA = MA./repmat(sqrt(sum(MA.^2,1)),ncha*lwin,1);
MB = MB./repmat(sqrt(sum(MB.^2,1)),ncha*lwin,1);


if strcmp(phase,'P')
    
    % -- First get CCs, flip polarity of wvfB where CC < 0
    CCF = sum(MA.*MB,1);
    MB  = MB.*repmat(sign(CCF),ncha*lwin,1);

    % -- Get 1-norm misfit, convert to number comparable to CC
    N1 = sum(abs(MA-MB),1); 
    N1 = 1-N1/max(N1);
     
    % -- For every peak, find indices of adjacent zero-crossings
    % -- This assumes the CCF will cross zero between every peak...
    sgnccf = sign(CCF);
    dsgn   = [0,diff(sgnccf)];
    zrcr   = [1,find(dsgn),nlg+1];
elseif strcmp(phase,'S')

    % -- Initialize polarity multiplier
    pm = cell(ncha,1);

    % -- Get sign of each CC (for each channel)
    for icha = 1:ncha
        jjc      = lwin*(icha-1)+[1:lwin];
        pol      = sign(sum(MA(jjc,:).*MB(jjc,:),1)) ;
        pm{icha} = repmat(pol,lwin,1);
    end
    
    pm  = vertcat(pm{:});
    CCF = sum(MA.*MB.*pm,1);
    N1  = sum(abs(MA-MB.*pm),1);
    N1  = 1-N1/max(N1);

    % -- For every peak, find indices of adjacent 0.5-crossings
    % -- This is altered from P-wave version for pure-positve CCF
    sgnccf = sign(CCF-0.5);
    dsgn   = [0,diff(sgnccf)];
    zrcr   = [1,find(dsgn),nlg+1];
    
else
    % -- 2020-05-12
    % -- ADDED THIS "GENERAL CASE FOR SIMPLE CC (POSITIVE ONLY)
    % -- Not 100% sure the "finding all the peaks stuff" below works
    
    % -- First get CCs, flip polarity of wvfB where CC < 0
    CCF = sum(MA.*MB,1);
    %MB  = MB.*repmat(sign(CCF),ncha*lwin,1);

    % -- Get 1-norm misfit, convert to number comparable to CC
    N1 = sum(abs(MA-MB),1); 
    N1 = 1-N1/max(N1);
     
    % -- For every peak, find indices of adjacent zero-crossings
    % -- This assumes the CCF will cross zero between every peak...
    sgnccf = sign(CCF);
    dsgn   = [0,diff(sgnccf)];
    zrcr   = [1,find(dsgn),nlg+1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Find minima of L1 misfit

[n1mx,icid] = max(abs(N1(2:end-1)));
cid         = icid+1;
[DT,CC]     = max_quad_interp(lg([cid-1:cid+1])*dt,CCF(cid-1:cid+1));

% -- List peaks
pks = zeros(nmax,2);

% -- Perhaps I could subtract a median from the CCF first?
pks(1,:) = [DT,CC];
j1 = zrcr(find(zrcr < cid,1,'last'));
j2 = zrcr(find(zrcr > cid,1))-1;
N1(j1:j2) = 0;

for ii = 2:nmax

    % -- Get next highest peak
    [n1mx,icid] = max(abs(N1(2:end-1)));
    cid         = icid+1;
    [pks(ii,1),pks(ii,2)] = max_quad_interp(lg([cid-1:cid+1])*dt,CCF(cid-1:cid+1));

    % -- Zero out this peak
    j1 = zrcr(find(zrcr < cid,1,'last'));
    j2 = zrcr(find(zrcr > cid,1))-1;
    N1(j1:j2) = 0;

    if ~nnz(N1)
        return
    end

end



