function DCplotter(sdr,x0,y0,r0,color,lw,hemi)
% function DCplotter(sdr,x,y,r0,color,lw,hemi)
%
% 2020-06-05
% This function produces a DC focal mechanism beachball plot 
% (assumes lower hemisphere) at a given location/size on an exiting plot.
% Can take in an arbitrary number of events, but the number of x0/y0/r0
% must match the numberof sdr
%
%     INPUTS
%
%      sdr = [strike1 dip1 rake1 strike2 dip2 rake2]
%       x0 = x-position on plot
%       y0 = y-position on plot
%       r0 = radius of beachball on plot
%    color = fill color of compressive quadrants e.g. 'k' or [0.5 0.5 0.5]
%       lw = linewidth
%     hemi = 'lower' (default) or 'upper' hemisphere
%
%   OUTPUTS
%
%      


hold on

Ne  = size(sdr,1);

% -- If we wanted to use upper-hemisphere projection, 
% -- this would have be -1
if nargin < 7
    hemi = 'lower';
    if nargin < 6
        lw = 2;
        if nargin < 5
            color = 'b';
            if nargin < 4
                x0 = 0;
                y0 = 0;
                r0 = 1;
                if Ne > 1
                    disp('Error: Multiple events input without x,y, and r!')
                    return
                end
            end
        end
    end
end
    
if strcmp(hemi,'upper')
	pm = -1;
else
    pm = 1;
end


%rng = [1:-0.01:-1]';
rng = [1:-0.005:0.98,0.97:-0.01:0.9,0.88:-0.02:0.02];
rng = [rng,0,-fliplr(rng)]';
phic = [0:2:360]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Outer circle
xc = cosd(phic);
yc = sind(phic);

for ii = 1:Ne

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- First nodal plane
    s = sdr(ii,1);
    d = sdr(ii,2);

    sv = [cosd(90-s),sind(90-s),0];
    dv = [cosd(-s)*cosd(d),sind(-s)*cosd(d),-sind(d)];

    ws = repmat(rng,1,3);
    wd = sqrt(1-ws.^2);
    V = ws.*sv + wd.*dv;

    x1 = pm*sqrt(1./(1-V(:,3))).*V(:,1);
    y1 = pm*sqrt(1./(1-V(:,3))).*V(:,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Second nodal plane
    s = sdr(ii,4);
    d = sdr(ii,5);

    sv = [cosd(90-s),sind(90-s),0];
    dv = [cosd(-s)*cosd(d),sind(-s)*cosd(d),-sind(d)];

    ws = repmat(rng,1,3);
    wd = sqrt(1-ws.^2);
    V  = ws.*sv + wd.*dv;

    x2 = pm*sqrt(1./(1-V(:,3))).*V(:,1);
    y2 = pm*sqrt(1./(1-V(:,3))).*V(:,2);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Fill compressive regions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Pure Reverse fault
    if sdr(ii,3) == 90
    
        xf1 = [x1; x2];
        yf1 = [y1; y2];
        xf2 = [];
        yf2 = [];
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Pure Normal fault    
    elseif sdr(ii,3) == -90
       
        phir = atan2d(y1(end),x1(end))+[2:2:178]';
        xf1 = [x1; cosd(phir)];
        yf1 = [y1; sind(phir)];
    
        phir = atan2d(y2(end),x2(end))+[2:2:178]';
        xf2 = [x2; cosd(phir)];
        yf2 = [y2; sind(phir)];
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- General case (rakes ~= +/- 90)    
    else    
    
        no1 = [cosd(-sdr(ii,1)).*sind(sdr(ii,2)), sind(-sdr(ii,1)).*sind(sdr(ii,2)), cosd(sdr(ii,2))];
        no2 = [cosd(-sdr(ii,4)).*sind(sdr(ii,5)), sind(-sdr(ii,4)).*sind(sdr(ii,5)), cosd(sdr(ii,5))];
        ax3 = cross(no1,no2);

        ax3 = -sign(ax3(3))*ax3/norm(ax3);

        xI = pm*sqrt(1/(1-ax3(3)))*ax3(1);
        yI = pm*sqrt(1/(1-ax3(3)))*ax3(2);

        [~,jm1] = min( (x1-xI).^2 + (y1-yI).^2 );
        [~,jm2] = min( (x2-xI).^2 + (y2-yI).^2 );
        
        if (jm1>1)*(jm1<length(x2))      
            if (x1(jm1+1)-xI)^2 + (y1(jm1+1)-yI)^2 < (x1(jm1-1)-xI)^2 + (y1(jm1-1)-yI)^2 
                x1 = [x1(1:jm1); xI; x1(jm1+1:end)];
                y1 = [y1(1:jm1); yI; y1(jm1+1:end)];
                jm1 = jm1+1;
            else
                x1 = [x1(1:jm1-1); xI; x1(jm1:end)];
                y1 = [y1(1:jm1-1); yI; y1(jm1:end)];
            end
        elseif jm1 > 1
            x1 = [x1; xI];
            y1 = [y1; yI];
            jm1 = jm1+1;
        else
            x1 = [xI; x1];
            y1 = [yI; y1];
        end
        
        if (jm2>1)*(jm2<length(x2))  
            if (x2(jm2+1)-xI)^2 + (y2(jm2+1)-yI)^2 < (x2(jm2-1)-xI)^2 + (y2(jm2-1)-yI)^2
                x2 = [x2(1:jm2); xI; x2(jm2+1:end)];
                y2 = [y2(1:jm2); yI; y2(jm2+1:end)];
                jm2 = jm2+1;
            else
                x2 = [x2(1:jm2-1); xI; x2(jm2:end)];
                y2 = [y2(1:jm2-1); yI; y2(jm2:end)];
            end
        elseif jm2 > 1
            x2 = [x2; xI];
            y2 = [y2; yI];
            jm2 = jm2+1;
        else
            x2 = [xI; x2];
            y2 = [yI; y2];
        end
        
        % -- How to identify compressional quadrants??
        % -- For now just assume something...
        if sdr(ii,3) >= 0
            xf1 = [x1(1:jm1); x2(jm2+1:end)];
            yf1 = [y1(1:jm1); y2(jm2+1:end)];
            xf2 = [x2(1:jm2); x1(jm1+1:end)];
            yf2 = [y2(1:jm2); y1(jm1+1:end)];
        else
            xf1 = [x1(1:jm1); x2(jm2-1:-1:1)];
            yf1 = [y1(1:jm1); y2(jm2-1:-1:1)];
            xf2 = [x2(end:-1:jm2); x1(jm1+1:end)];
            yf2 = [y2(end:-1:jm2); y1(jm1+1:end)];
        end

        phiE = atan2d(yf1([end,1]),xf1([end,1]));
        if (phiE(1)<phiE(2))*(phiE(2)-phiE(1)<=180)
            phir = [phiE(1)+2:2:phiE(2)]';
        elseif (phiE(1)<phiE(2))*(phiE(2)-phiE(1)>180)
            phir = [phiE(1)-2:-2:phiE(2)-360]';
        elseif (phiE(1)>=phiE(2))*(phiE(1)-phiE(2)<=180)
            phir = [phiE(1)-2:-2:phiE(2)]';
        elseif  (phiE(1)>=phiE(2))*(phiE(1)-phiE(2)>180)
            phir = [phiE(1)+2:2:phiE(2)+360]';
        end
        
        xf1 = [xf1; cosd(phir)];
        yf1 = [yf1; sind(phir)];
       

        phiE = atan2d(yf2([end,1]),xf2([end,1]));
        if (phiE(1)<phiE(2))*(phiE(2)-phiE(1)<=180)
            phir = [phiE(1)+2:2:phiE(2)]';
        elseif (phiE(1)<phiE(2))*(phiE(2)-phiE(1)>180)
            phir = [phiE(1)-2:-2:phiE(2)-360]';
        elseif (phiE(1)>=phiE(2))*(phiE(1)-phiE(2)<=180)
            phir = [phiE(1)-2:-2:phiE(2)]';
        elseif  (phiE(1)>=phiE(2))*(phiE(1)-phiE(2)>180)
            phir = [phiE(1)+2:2:phiE(2)+360]';
        end
        xf2 = [xf2; cosd(phir)];
        yf2 = [yf2; sind(phir)];
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Actual plot commands
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -- Should add "screen" or "print" options...
    
    fill(r0(ii)*xc +x0(ii), r0(ii)*yc +y0(ii), 'w')
    plot(r0(ii)*xc +x0(ii), r0(ii)*yc +y0(ii), 'k','Linewidth',lw);
    fill(r0(ii)*xf1+x0(ii), r0(ii)*yf1+y0(ii), color)
    fill(r0(ii)*xf2+x0(ii), r0(ii)*yf2+y0(ii), color)
    
    
    
    %plot(r0(ii)*x1 +x0(ii), r0(ii)*y1 +y0(ii), 'k','LineWidth',lw)
    %plot(r0(ii)*x2 +x0(ii), r0(ii)*y2 +y0(ii), 'k','LineWidth',lw)
    %plot(r0(ii)*xc +x0(ii), r0(ii)*yc +y0(ii), 'k','Linewidth',lw);
end



    %jp = find( (phic>min(phiENDS)).*(phic<max(phiENDS)) );
    %jp = [ find( (phic>min(phiENDS)).*(phic<max(phiENDS)) );
    
    %{
    if phiENDS(1) < phiENDS(2)
        phir = [phiENDS(1)+2:2:phiENDS(2)+360];
        xf1 = [xf1; r0(ii)*cosd(phir)+x0(ii)];
        yf1 = [yf1; r0(ii)*cosd(phir)+x0(ii)];
    else
        xf1 = [xf1; flipud(xc(jp))];
        yf1 = [yf1; flipud(yc(jp))];
    end
    
    jp = find( (phic>min(phiENDS)).*(phic<max(phiENDS)) );
    if phiENDS(1) < phiENDS(2)
        xf2 = [xf2; xc(jp)];
        yf2 = [yf2; yc(jp)];
    else
        xf2 = [xf2; flipud(xc(jp))];
        yf2 = [yf2; flipud(yc(jp))];
    end
    
    %}


