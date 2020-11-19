function plotFirstMotions(mp,mn,hemi,ue,fm,phaprms)
% function plotFirstMotions(mp,mn,hemi,ue,fm,phaprms)
%
% 2018-09-12
% This function adds first-motions to MT plots.
%
%
%  INPUTS
%
%         mp == marker colour, style for positive (up) first motions      
%         mn == marker colour, style for negative (down) first motions 
%         hemi == 'lower' (default) or 'upper' hemisphere 
%                  NOT ACTUALLY DOING ANYTHING RIGHT NOW
%         ue == event numbers being plotted
%         fm == matrix containing first-motion picks [EV,STA,POL,PHA]
%   pharprms == matrix containing phase azimuth, take-off angles, distances
%
%

ne  = length(ue);

% -- First need to organize subplots.
npx = fix(ceil(sqrt(ne)));
npy = fix(ceil(ne/npx));

%[axs,pos] = tightsubplot(npx,npy,0.01,0.01,0.01);
%fig  = gcf;

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
   
    axes(axs(ie));    
    hold on

   
    %fmE = fm(fm(:,1)==ue(ie),:);
    %fmP = fmE(fmE(:,3)==1,:);
    %fmN = fmE(fmE(:,3)==-1,:);
    %phaP = fmP(:,4);
    %phaN = fmN(:,4);
    fmE = fm(fm(:,1)==ue(ie),:);
    fmP = fmE(fmE(:,4)==1,:);
    fmN = fmE(fmE(:,4)==-1,:);
    phaP = fmP(:,3);
    phaN = fmN(:,3);
    
    % -- Get azimuths and take-off angles
    aziP = phaprms(phaP,1);
    aziN = phaprms(phaN,1);
    %if strcmp(hemi,'lower')
    %    aziP = aziP+180;
    %    aziN = aziN+180;
    %    aziP(aziP>360) = aziP(aziP>360)-360;
    %    aziN(aziN>360) = aziN(aziN>360)-360;
    %end
    toaP = phaprms(phaP,2);
    toaN = phaprms(phaN,2);

    ig90P = find(toaP > 90);
    il90P = find(toaP < 90);
    ig90N = find(toaN > 90);
    il90N = find(toaN < 90);
    
    % -- Plot distances
    rP        = zeros(length(aziP),1);
    rN        = zeros(length(aziN),1);
    rP(ig90P) = sqrt(2)*cosd(toaP(ig90P)/2);
    rN(ig90N) = sqrt(2)*cosd(toaN(ig90N)/2);
    rP(il90P) = sqrt(2)*cosd((180-toaP(il90P))/2);
    rN(il90N) = sqrt(2)*cosd((180-toaN(il90N))/2);

    % -- Convert azimuths, plot distances to plot x,y
    xP = rP.*cosd(90-aziP);
    yP = rP.*sind(90-aziP);
    xN = rN.*cosd(90-aziN);
    yN = rN.*sind(90-aziN);

    keyboard
    % -- Plot phases 
    if strcmp(hemi,'lower')
        xP(il90P) = -xP(il90P);
        yP(il90P) = -yP(il90P);
        xN(il90N) = -xN(il90P);
        yN(il90N) = -yN(il90P);
    else
        xP(ig90P) = -xP(ig90P);
        yP(ig90P) = -yP(ig90P);
        xN(ig90P) = -xN(ig90P);
        yN(ig90P) = -yN(ig90P);
    end
    plot(xP(il90P),yP(il90P),mp)
    plot(xN(il90N),yN(il90N),mn)
    %plot(xsP(ig90P),ysP(ig90P),'bo')
    %plot(xsS(ig90S),ysS(ig90S),'rv')
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
end
