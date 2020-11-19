function  MTspherePS(M,clr,quiver)
% function  MspherePS(M,clr,quiver)
%
% 2018-09-07
% This function makes a 3D sphere plot of a Moment tensor. The default behavior
% is to show the P field with a red-white-blue colormap. And show the S field
% with arrows.
%
%  INPUTS
%
%        M == a 3x3 moment tensor in ENU convention, or a 6x1 vector of
%             [M11,M22,M33,M12,M13,M23]
%      clr == 'b' for compressional region color. Input empty string '' for no
%              P field or 'rwb' for more detailed colormap showing 
%              amplitude intensity.
%   quiver == 1 (default) to include arrows for S field, 0 otherwise
%
%

if nargin < 2
    clr = 'rwb';
elseif clr == 'rwb'
    r1=[(0:127)/128,ones(1,128)];
    g1=[(0:127)/127,(127:-1:0)/127];
    b1=[ones(1,128),(127:-1:0)/127];
    rwb=flipud([r1',g1',b1']);
    colormap(rwb)
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
    colormap([1 1 1; clr])
end

if nargin < 3
    quiver = 1;
end

if length(M) == 3
    M = M/sqrt(sum(M(:).^2)/2);
    M = M([1 5 9 2 3 6]');
else
    M = M/sqrt(sum(M(1:3).^2+2*M(4:6).^2)/2);
end 

% -- Define resolution (64, 128, or 256)
res = 256;

% -- Define points of sphere
[x,y,z] = sphere(res);

% -- Find azimuth, radius to each grid point. 
phi = atan2d(y,x);
r   = sqrt(x.^2 + y.^2);
tht = 90-atan2d(z,r);

% -- Use take off angles to compute coefficients for each of the 6 MT parts
gx = cosd(phi).*sind(tht);
gy = sind(phi).*sind(tht);
gz = cosd(tht);

gx = gx(:);
gy = gy(:);
gz = gz(:);
g = [gx.^2, gy.^2, gz.^2, 2*gx.*gy, 2*gx.*gz, 2*gy.*gz];

% -- Compute amplitudes at each angle
% -- Set default as something very small and negative so that part outside
% -- the circle plots as white in 'b' mode.
amp = g*M;
amp = reshape(amp,size(phi));

clf
hold on


if clr == 'rwb'
    surf(x,y,z,amp,'edgecolor','none')
elseif length(clr)
    
    
    surf(x,y,z,amp,'edgecolor','none')
   
end

caxis([-max(abs(amp(:))),max(abs(amp(:)))])

% Plot equators
r1 = res/4+1;
r2 = res/2+1;
r3 = res*(3/4)+1;
plot3(x(:,1),y(:,1),z(:,1),'k')
plot3(x(:,r1),y(:,r1),z(:,r1),'k')
plot3(x(:,r2),y(:,r2),z(:,r2),'k')
plot3(x(:,r3),y(:,r3),z(:,r3),'k')
plot3(x(r2,:),y(r2,:),z(r2,:),'k')

% If adding quiver plot, redo grid at coarser scale
if quiver

    [x,y,z] = sphere(20);

    phi = atan2d(y,x);
    r   = sqrt(x.^2 + y.^2);
    tht = 90-atan2d(z,r);

    gx = cosd(phi).*sind(tht);
    gy = sind(phi).*sind(tht);
    gz = cosd(tht);

    gx = gx(:);
    gy = gy(:);
    gz = gz(:);
    gP = [gx.^2, gy.^2, gz.^2, 2*gx.*gy, 2*gx.*gz, 2*gy.*gz];

    zr  = zeros(size(gx));
    gX = [gx, zr, zr, gy, gz, zr] - repmat(gx,1,6).*gP;
    gY = [zr, gy, zr, gx, zr, gz] - repmat(gy,1,6).*gP;
    gZ = [zr, zr, gz, zr, gx, gy] - repmat(gz,1,6).*gP;

    ampX = gX*M;
    ampY = gY*M;
    ampZ = gZ*M;
    ampX = reshape(ampX,size(phi));
    ampY = reshape(ampY,size(phi));
    ampZ = reshape(ampZ,size(phi));

    %surf(x,y,z,'edgecolor',[0.5 0.5 0.5])
    quiver3(1.03*x,1.03*y,1.03*z,ampX,ampY,ampZ,'k','LineWidth',2,'MaxHeadSize',0.8)
end

axis equal
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
xlim([-1.2,1.2])
ylim([-1.2,1.2])
zlim([-1.2,1.2])
xlabel('X')
ylabel('Y')
zlabel('Z')

