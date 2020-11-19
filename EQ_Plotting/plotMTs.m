function plotMTs(M,clrmd,hemi,filled)
% function plotMTs(M,clrmd,hemi,filled)
%
% 2018-05-15
% This file plots all moment tensors stored in the cell array 'M'.
% It also finds which station were used for that event, their take-off angles
% and azimuths, and add to plot.
% raP  = [EV1,EV2,STA,INS,RA,RM1,RM2,P1,P2];
% raS3 = [evA,evB,evC,STA,INS, (1:5)
%         Babc,Bacb,RM,CCabc   (6:9)
%         CCab,CCab,CCbc,      (10:12)
%         Pa,Pb,Pc]            (13:15)
%
% 2018-05-28 
% Updated to be a function
%
%   INPUTS
%
%        M  == a cell array of moment tensors in ENU convention      
% 
%     clrmd == 'rwb' for more detailed colormap showing amplitude intensity.
%              any color code for compressive regions to be that color--vector
%              graphics used in this case. 
%      hemi == 'lower' (default) or 'upper' hemisphere
%    filled == 1 (default) or 0 .Whether or not to fill compressive region 
%              (doesn't apply if using 'rwb' colormode). If 0, then 'clr' will 
%              apply to the contour.
%

%set(gcf,'Visible','Off');

if nargin < 4
    filled = 1;
end

ne  = length(M);
axs = findall(gcf,'type','axes');

% -- First need to organize subplots.
npx = fix(ceil(sqrt(ne)));
npy = fix(ceil(ne/npx));
npt = npx*npy;

if length(axs)
    
    % -- Get axis positions, sort
    xya = zeros(ne,2);
    for ii = 1:length(axs)
        xya(ii,:) = axs(ii).Position([2,1]);
    end
    
    xya(:,1) = 1-xya(:,1);
    [~,jj]   = sortrows(xya);
    axs = axs(jj);

    %if ne < npt
    %    axs(1:(npt-ne)) = [];
    %end
    %axs = flipud(axs);
else
    [axs,pos] = tightsubplot(npx,npy,0.01,0.01,0.01);
end

for ie = 1:ne

    % -- Allow only some subplots to be drawn on by checking if an MT was
    % -- provided for this slot
    if ~length(M{ie})
        continue
    end
    % -- Set up axes
    %keyboard
    axes(axs(ie));    
    %ax = axs(ie);
    hold on

    % -- Add MT plot
    %focalmech(M{ie},0,0,2);
    MTplotter(M{ie},clrmd,hemi,filled);
  
    axis equal
    xlim([-1.1,1.1])
    ylim([-1.1,1.1])
    ax.XTick = [];
    ax.YTick = [];
   
end
