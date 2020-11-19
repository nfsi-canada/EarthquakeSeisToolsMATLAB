function plotNodalPlane(sdr,linespec,hemi)
% function plotNodalPlane(sdr,linespec,hemi)
%
% 2020-06-04
% This function plots a nodal plane on a beach ball plot given a 
% strike/dip/rake.
%
%    INPUTS
%
%      sdr  = [strike,dip,rake] if more than one row is provided
%             each plane is plotted
%  linespec = 
%      hemi = 'lower' (default) or 'upper' hemisphere



if nargin < 3
    hemi = 'lower';
    if nargin < 2
        linespec = 'k';
    end    
end


if strcmp(hemi,'upper')
    pm = -1;
else
    pm = 1;
end


for ii = 1:size(sdr,1)
    s = sdr(ii,1);
    d = sdr(ii,2);

    sv = [cosd(90-s),sind(90-s),0];
    dv = [cosd(-s)*cosd(d),sind(-s)*cosd(d),-sind(d)];

    ws = repmat([-1:0.005:1]',1,3);
    wd = sqrt(1-ws.^2);
    V = ws.*sv + wd.*dv;

    x = pm*sqrt(1./(1-V(:,3))).*V(:,1);
    y = pm*sqrt(1./(1-V(:,3))).*V(:,2);
    
    plot(x,y,linespec,'LineWidth',2.5)
    
end


