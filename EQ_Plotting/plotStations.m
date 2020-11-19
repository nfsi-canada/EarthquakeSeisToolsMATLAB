function   plotStations(M,clrmd,hemi,ue,evsM,raP,raS3,phaprms)
% function plotStations(M,clrmd,hemi,ue,evsM,raP,raS3,phaprms)
%
% 2018-11-06
% Rough changes made from plotMTandStations to ignore MTs, plot angle contours
%
% 2018-05-28 
% Updated to be a function
%
%   INPUTS
%
%        M  == a cell array of moment tensors in ENU convention      
% 
%     clrmd == 'b' for typical black and white plot, or 'rwb' for more detailed
%               colormap showing amplitude intensity. Default is 'rwb';
%     hemi  == 'lower' (default) or 'upper' hemisphere
%      
%

evsM = find(ismember(ue,evsM));

ne  = length(evsM);

% -- First need to organize subplots.
npx = fix(ceil(sqrt(ne)));
npy = fix(ceil(ne/npx));

axs = findall(gcf,'type','axes');

if length(axs)
    axs  = flipud(axs);
    swAx = 0;
else
    [axs,pos] = tightsubplot(npx,npy,0.01,0.01,0.01);
    swAx = 1;
    aa = [0:15:90];
    ab = [0:30:150];
end

for ie = 1:ne

    je = evsM(ie);
    
    % -- Set up axes
    %ax = subplot(npx,npy,ie);
    axes(axs(ie));    
    hold on

    % -- Add MT plot
    %focalmech(M{ie},0,0,2);
    %MTplotter(M{je},clrmd,hemi)
    
    % -- Station stuff
    tmpAP = raP(raP(:,1)==je,[3, 9]);
    tmpBP = raP(raP(:,2)==je,[3,10]);
    tmpP  = unique([tmpAP; tmpBP],'rows');
    staP  = tmpP(:,1);
    phaP  = tmpP(:,2);
    ncmpP = hist([tmpAP(:,1); tmpBP(:,1)],staP);

    tmpAS = raS3(raS3(:,1)==je,[4,13]);
    tmpBS = raS3(raS3(:,2)==je==je,[4,14]);
    tmpCS = raS3(raS3(:,3)==je,[4,15]);
    tmpS  = unique([tmpAS; tmpBS; tmpCS],'rows');
    staS  = tmpS(:,1);
    phaS  = tmpS(:,2);
    ncmpS = hist([tmpAS(:,1); tmpBS(:,1); tmpCS(:,1)],staS);

    %staC = intersect(staP,staS);
    
    % -- Get azimuths and take-off angles
    aziP = phaprms(phaP,1);
    aziS = phaprms(phaS,1);
    toaP = phaprms(phaP,2);
    toaS = phaprms(phaS,2);

    % -- 
    ig90P = find(toaP > 90);
    il90P = find(toaP < 90);
    ig90S = find(toaS > 90);
    il90S = find(toaS < 90);
    
    % -- Plot distances
    rP        = zeros(length(staP),1);
    rS        = zeros(length(staS),1);
    rP(ig90P) = sqrt(2)*cosd(toaP(ig90P)/2);
    rS(ig90S) = sqrt(2)*cosd(toaS(ig90S)/2);
    rP(il90P) = sqrt(2)*cosd((180-toaP(il90P))/2);
    rS(il90S) = sqrt(2)*cosd((180-toaS(il90S))/2);

    % -- Convert azimuths, plot distances to plot x,y
    xsP = rP.*cosd(90-aziP);
    ysP = rP.*sind(90-aziP);
    xsS = rS.*cosd(90-aziS);
    ysS = rS.*sind(90-aziS);

    % -- Now actually plot stations
    %scatter(xsP(il90P),ysP(il90P),38,ncmpP(il90P),'o','filled')
    %scatter(xsS(il90S),ysS(il90S),38,ncmpS(il90S),'s','filled')
    %scatter(xsP(ig90P),ysP(ig90P),38,ncmpP(ig90P),'o')
    %scatter(xsS(ig90S),ysS(ig90S),38,ncmpS(ig90S),'s')
   
    % -- Plot contours on angle
    if swAx
        for ii = 1:length(aa)
            toa = aa(ii);
            azi = 1:360;
            r   = sqrt(2)*cosd((180-toa)/2);
            x   = r.*cosd(90-azi);
            y   = r.*sind(90-azi);
            plot(x,y,'Color',[0.7,0.7,0.7])
        end
        for ii = 1:length(ab)
            toa = -90:90;
            azi = ab(ii);
            r   = sqrt(2)*cosd((180-toa)/2);
            x   = r.*cosd(90-azi);
            y   = r.*sind(90-azi);
            plot(x,y,'Color',[0.7,0.7,0.7])
        end
    end

    plot(xsP(il90P),ysP(il90P),'g+','MarkerSize',2)
    plot(xsS(il90S),ysS(il90S),'rv','MarkerSize',2,'Linewidth',0.1)
    
    %plot(xsP(il90P),ysP(il90P),'go','MarkerFaceColor','g','MarkerSize',2)
    %plot(xsS(il90S),ysS(il90S),'rv','MarkerFaceColor','r','MarkerSize',2)
    %plot(xsP(ig90P),ysP(ig90P),'go','MarkerSize',2)
    %plot(xsS(ig90S),ysS(ig90S),'rv','MarkerSize',2)

    axis equal
    xlim([-1.1,1.1])
    ylim([-1.1,1.1])
    ax.XTick = [];
    ax.YTick = [];
end
