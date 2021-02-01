function  A = RayPath2NodeWeights3D(gx,gy,gz,RP)
% function A = RayPath2NodeWeights3D(gx,gy,gz,RP)
%
% 2020-06-14
% Takes in discretized ray-paths and to computes distances travelled 
% through each cell... really it's best to think of the velocity model
% as nodes. 
%
%    INPUTS
% 
%    gx,gy,gz = vectors of model node coordinates
%          RP = a cell-array containing Nray lists of ray-path
%               coordiantes [x,y,z].
%
%   OUTPUTS
%
%           A = Nray x Ncell matrix indicating distance per cell
%
%          %%%D = a cell array containing lists of node distance weights
%          %%%    [node#, distance], where node# 


Nray = length(RP);
%D = cell(Nray,1);

dgx = diff(gx(1:2))/2;
dgy = diff(gy(1:2))/2;
dgz = diff(gz(1:2))/2;

%mindg = min([dgx,dgy,dgz]);
%[GX,GY,GZ] = ndgrid(gx,gy,gz);


% -- There must be a better way, but the is the best/fastest way I
% -- know how to get linear indices from ranges in this direction involves
% -- This extra 3D matrix
Nx = length(gx);
Ny = length(gy);
Nz = length(gz);
Nc = Nx*Ny*Nz;
Mj = repmat((1:Nx)',1,Ny,Nz) + repmat(Nx*(0:Ny-1),Nx,1,Nz)  ...
     + repmat( permute((Nx*Ny)*(0:Nz-1),[3 1 2]),Nx,Ny,1);

A = sparse(Nray,Nc);

% -- Run in parallel if pool is open
pool = gcp('nocreate');
if isempty(pool)
    Nworker = 0;
else
    Nworker = pool.NumWorkers;
end

col = @(x) x(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor (iray = 1:Nray, Nworker);
    
    % -- What if ray-path steps are far apart??
    X = RP{iray};
    
    dseg = vecnorm(diff(X)')';
    Nseg = length(dseg);
    
    d_node = sparse(Nc,1);
    
    for ii = 1:Nseg
        
        Xb = sort(X(ii:ii+1,:));

        % -- 
        jx = find((gx > Xb(1,1)-dgx).*(gx < Xb(2,1)+dgx));
        jy = find((gy > Xb(1,2)-dgy).*(gy < Xb(2,2)+dgy));
        jz = find((gz > Xb(1,3)-dgz).*(gz < Xb(2,3)+dgz));
        
        jn = col(Mj(jx,jy,jz));

        % -- Convert 3D grid indices to indices of gx,gy,gz vectors
        % -- This is slightly slower than referencing GX,GY,GZ grids, but
        % -- I think it's better not to have to pass large grids to workers
        jz = ceil(jn/(Nx*Ny));
        jy = ceil((jn-Nx*Ny*(jz-1))/Nx);
        jx = jn - (jy-1)*Nx - (jz-1)*Nx*Ny;

        dn2 = sum([gx(jx)-Xb(1,1),gy(jy)-Xb(1,2),gz(jz)-Xb(1,3), ... 
                   gx(jx)-Xb(2,1),gy(jy)-Xb(2,2),gz(jz)-Xb(2,3)].^2,2)
        
   
        wn = 1./dn2;
        
        d_node(jn) = d_node(jn) + dseg(ii)*wn/sum(wn);        
    
    end
    
    
    A(iray,:) = d_node;
    %jn = find(d_node);
    %A(iray,jn) = d_node(jn);
    %D{iray} = [jn,d_node(jn)];
    
    
end

%jx = find((gx > Xb(1,1)-dgx).*(gx < Xb(2,1)+dgx));
% -- I think this find statement is too slow...
%jn = find( (GX > Xb(1,1)-dgx).*(GX < Xb(2,1)+dgx).* ...
%           (GY > Xb(1,2)-dgy).*(GY < Xb(2,2)+dgy).* ...
%           (GZ > Xb(1,3)-dgz).*(GZ < Xb(2,3)+dgz) );

% -- Weights must sum to 1...relate to mean distance from ray-path
% -- end poitns??
%dn2 = sum([GX(jn)-Xb(1,1),GY(jn)-Xb(1,2),GZ(jn)-Xb(1,3), ... 
%           GX(jn)-Xb(2,1),GY(jn)-Xb(2,2),GZ(jn)-Xb(2,3)].^2,2); 

