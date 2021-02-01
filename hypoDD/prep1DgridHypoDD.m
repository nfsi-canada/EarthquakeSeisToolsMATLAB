function [Ft,Fg,Fv] = prep1DgridHypoDD(model,Xmax,Zmax)

% 2021-01-14
% This function prepares matrices of time, take-off angle, and source
% velocity to each station (if they different elevations). If station
% elevations aren't used then gT, gA can be the same for all stations...
%
% Is it better to assume first layer extends up to elevation???
%    Or shift all layers upward?? ... let's do the shift all layers...
%    then the model is always the same
%
% INPUTS
%
%   model = [Z(km), V(km s-1)] N_layer x 2 matrix
%    Xmax = maximum epicentral distance
%    Zmax = maximum event depth
%
% OUTPUTS
%
%      Ft = interpolant of travel time given distance,depth
%      Fg = interpolant of vertical direction cosine (positive=down)
%           given distance,depth
%      Fv = interpolant of source velocity given depth


x = exp(linspace(0.1,log(Xmax),100))-1;
z = exp(linspace(0.1,log(Zmax),100))-1;    
Nx = length(x);
Nz = length(z);

% -- Make 2D grid
[X,Z] = ndgrid(x,z);

t = zeros(Nx,Nz);
a = zeros(Nx,Nz);

% -- Do ray tracing
for iz = 1:Nz
    [t(:,iz),a(:,iz)] = RayTrace(z(iz),x,model);
end

% -- Convert take-off angles to vertical direction cosine
% -- (with positive == DOWN)
gz = -cosd(a);

% -- Make interpolants for travel-time and direction cosines
Ft = scatteredInterpolant(X(:),Z(:),t(:));
Fg = scatteredInterpolant(X(:),Z(:),gz(:));

% -- Make interpolant from source velocities
%Fv = scatteredInterpolant(model(:,1),model(:,2));
Fv = @(z) interp1(model(:,1),model(:,2),z);