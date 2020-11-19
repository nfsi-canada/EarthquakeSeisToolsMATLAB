function AddFirstMotions(FM,mp,mn,x0,y0,r0,lw,hemi)
% function AddFirstMotions(mp,mn,hemi,ue,fm,phaprms)
%
% 2018-09-12
% This function adds first-motions to an MT plot.
%
%
%  INPUTS
%
%         FM == matrix containing [azimuth,take-off angle, and
%                                  +/- 1 indicating first motion]                            
%         mp == marker colour, style for positive (up) first motions      
%         mn == marker colour, style for negative (down) first motions 
%         x0 == center of focal mechanism x-coordinate in plot (default 0)
%         y0 == center of focal mechanism y-coordinate in plot (default 0)
%         r0 == radius of focal mechanism circle in plot (default 1)
%         lw == linewidth (default == 1)
%       hemi == 'lower' (default) or 'upper' hemisphere 
%

if nargin < 8
    hemi = 'lower';
    if nargin < 7
        lw = 3;
        if nargin < 6
            x0 = 0;
            y0 = 0;
            r0 = 1;  
            if nargin < 3
                mp = 'g+';
                mn = 'ro';
            end
        end
    end
end
   
hold on

Nfm = size(FM,1);
azi = FM(:,1);
toa = FM(:,2);
jg90 = find(toa >= 90);
jl90 = find(toa < 90);   


% -- Plot distances
r       = zeros(Nfm,1);
r(jg90) = sqrt(2)*cosd(toa(jg90)/2);
r(jl90) = sqrt(2)*cosd((180-toa(jl90))/2);

% -- Convert azimuths, plot distances to plot x,y
x = r.*cosd(90-azi);
y = r.*sind(90-azi);
   
% -- Flip direction of phases leaving opposite hemisphere 
if strcmp(hemi,'lower')
    x(jl90) = -x(jl90);
    y(jl90) = -y(jl90);
else
    x(jg90) = -x(jg90);
    y(jg90) = -y(jg90);
end

x = x0 + x*r0;
y = y0 + y*r0;

% -- Get indices of positive and negative motions and plot
jP = find(FM(:,3) > 0);
jN = find(FM(:,3) < 0);
plot(x(jP),y(jP),mp,'LineWidth',lw)
plot(x(jN),y(jN),mn,'LineWidth',lw)
    