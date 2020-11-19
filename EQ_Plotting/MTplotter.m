function bb = MTplotter(M,clr,hemi,filled)
% function bb = MTplotter(M,clr,hemi,filled)
%
% This function makes a beach-ball plot of a Moment tensor. Blue are compressive
% (up-first) regions and red are tensional (down-first) regions. 
%
% 2018-07-03
% Updated to use vecotr graphics if not using 'rwb' mode. Also allows 
% compressive regions to be filled with any color.
%
%  INPUTS
%
%        M == a 3x3 moment tensor in ENU convention, or a 6x1 vector of
%             [M11,M22,M33,M12,M13,M23]
%      clr == 'rwb' for more detailed colormap showing amplitude intensity.
%              any color code for compressive regions to be that color--vector
%              graphics used in this case. 
%     hemi == 'lower' (default) or 'upper' hemisphere
%   filled == 1 (default) or 0 .Whether or not to fill compressive region 
%             (doesn't apply if using 'rwb' colormode). If 0, then 'clr' will 
%             apply to the contour.
%
%  OUTPUTS
%
%      bb == the plot.
%
%

sw_rwb = 0;

if nargin < 2
    clr = [0 0 0];
elseif clr == 'rwb'
    ;
elseif strcmp(class(clr),'char')
    cc = {'k', [0 0 0], 
    'b', [0 0 1],
    'c', [0 1 1], 
    'g', [0 1 0], 
    'm', [1 0 1], 
    'r', [1 0 0], 
    'y', [1 1 0]};

    jc  = cellstrfind(cc(:,1),clr);
    clr = cc{jc,2};
end

if nargin < 3
    hemi = 'lower';
end
if nargin < 4;
    filled = 1;
end

if length(M) == 3
    M = M/sqrt(sum(M(:).^2)/2);
    M = M([1 5 9 2 3 6]');
else
    M = M/sqrt(sum(M(1:3).^2+2*M(4:6).^2)/2);
end 

if clr == 'rwb'
    r1=[(0:127)/128,ones(1,128)];
    g1=[(0:127)/127,(127:-1:0)/127];
    b1=[ones(1,128),(127:-1:0)/127];
    rwb=flipud([r1',g1',b1']);
    colormap(rwb)
else
    colormap([1 1 1; clr])
end

% -- Should find pixel locations first, then find angles to them
% -- Just use square grid
%[px,py] = meshgrid([-1:0.005:1]);
% -- USE LOWER RES UNTIL PUBLISHING!!
[px,py] = meshgrid([-1:0.02:1]);

% -- Find azimuth, radius to each grid point. 
phi = atan2d(py,px);
r   = sqrt(px.^2 + py.^2);
jin = find(r < 1);
jpr = intersect(jin,find(r>0.996));

% -- Take-off angles (theta) depends on if plotting lower or upper hemisphere
if hemi == 'lower'
    tht = 2*acosd(r/sqrt(2));
elseif hemi == 'upper'
    tht = 360-2*acosd(r/sqrt(2));
end

% -- Use take off angles to compute coefficients for each of the 6 MT parts
gx = cosd(phi(jin)).*sind(tht(jin));
gy = sind(phi(jin)).*sind(tht(jin));
gz = cosd(tht(jin));

g = [gx.^2, gy.^2, gz.^2, 2*gx.*gy, 2*gx.*gz, 2*gy.*gz];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- This checks amplitudes around the perimeter to see if all are positive
% -- If so then the default is dark, with white on top. It's not perfect, if 
% -- there is a tensional zone that gets very close to the edge, then this might
% -- fail.
tmpamp = zeros(size(r));
tmpamp(jin) = g*M;
if ~length(find(tmpamp(jpr) > 0))
   amp   = -1e-6*ones(size(r));
   cmode = 0;
else
    amp   = 1e-6*ones(size(r));
    cmode = 1;
end
amp(jin) = g*M;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Plot MT, raster if 'rwb' mode, vector otherwise
hold on
if clr == 'rwb'
    imagesc([-1,1],[-1,1],amp);
else
    cntr = contour(px,py,amp,[0 0],'k-','Linewidth',1,'visible','off');
    
    if length(cntr)

        % -- Find number of compressive/tensional regions
        jcr = [1];
        ncr = 1;
        while 1
            jcr(end+1) = jcr(end)+cntr(2,jcr(end))+1;
            if jcr(end) > size(cntr,2)
                break
            end
            ncr = ncr+1;
        end

        % -- Different plot modes if compressional or extensional event
        if cmode == 1

            if filled
                % -- Fill circle (default is compressional)
                phic = 0:0.01:2*pi;
                xc   = cos(phic);
                yc   = sin(phic);
                fill(xc,yc,clr);

                for ii = 1:ncr
                    fill(cntr(1,jcr(ii)+1:jcr(ii+1)-1),cntr(2,jcr(ii)+1:jcr(ii+1)-1),'w');
                    plot(cntr(1,jcr(ii)+1:jcr(ii+1)-1),cntr(2,jcr(ii)+1:jcr(ii+1)-1),'k-','Linewidth',1);
                end
            else
                for ii = 1:ncr
                    plot(cntr(1,jcr(ii)+1:jcr(ii+1)-1),cntr(2,jcr(ii)+1:jcr(ii+1)-1),'Color',clr,'Linewidth',1);
                end
            end


        elseif cmode == 0 

            % -- If filled==1, fill compressive regions with vector polygon
            % -- then plot outline seperately, otherwise just plot outline
            if filled
                % -- Fill circle (default is tensional)
                phic = 0:0.01:2*pi;
                xc   = cos(phic);
                yc   = sin(phic);
                fill(xc,yc,'w')

                for ii = 1:ncr
                    fill(cntr(1,jcr(ii)+1:jcr(ii+1)-1),cntr(2,jcr(ii)+1:jcr(ii+1)-1),clr);
                    plot(cntr(1,jcr(ii)+1:jcr(ii+1)-1),cntr(2,jcr(ii)+1:jcr(ii+1)-1),'k-','Linewidth',1);
                end
            else
                for ii = 1:ncr
                    plot(cntr(1,jcr(ii)+1:jcr(ii+1)-1),cntr(2,jcr(ii)+1:jcr(ii+1)-1),'Color',clr,'Linewidth',1);
                end
            end
        end
    
    end
end

axis equal
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'visible','off')
xlim([-1,1])
ylim([-1,1])
box off

caxis([-max(abs(amp(:))),max(abs(amp(:)))])

% -- Add black circle (if filling in compressional portions, or using rwb)
if filled
    phic = 0:0.01:2*pi;
    xc   = cos(phic);
    yc   = sin(phic);
    plot(xc,yc,'k','Linewidth',1.5');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Compute amplitudes at each angle
% -- Set default as something very small and negative so that part outside
% -- the circle plots as white in 'b' mode.

% -- OR, if there is a negative isotropic component, plot whole thing as clr,
% -- then plot white on top (not in 'rwb' mode)
%{
ISO = sum(M(1:3));
if ISO <= 0 || strcmp(clr,'rwb')
    amp   = -1e-6*ones(size(r));
    cmode = 0;
else
    amp   = 1e-6*ones(size(r));
    cmode = 1;
end
amp(jin) = g*M;
%}
% -- Might have to check if MT is more 
%{
tmpamp = g*M;
if median(tmpamp) <= 0 
   amp   = -1e-6*ones(size(r));
   cmode = 0;
else
    amp   = 1e-6*ones(size(r));
    cmode = 1;
end
amp(jin) = g*M;
%}


