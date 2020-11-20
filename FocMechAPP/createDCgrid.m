% 2020-06-26
% This creates a ~uniform grid of DC moment tensors to quickly perform
% grid-searches later. The code actually creates the grid via random 
% rotations, so it is not perfect, but for N >~ 1e4 it's very close to 
% uniform

clear all

% -- How many focal mechanisms in grid?
N = 4e5;
M = uniform_grid_MT(N);

ng = length(M);

% -- Initialize matrices to store MTs, Str-Dip-Rke-Slp, plane normals
crsX = zeros(6,ng);
sdr  = zeros(ng,6);
xyzP = zeros(3,ng);
xyzT = zeros(3,ng);
xyzI = zeros(3,ng);
indx = [1 5 9 2 3 6]';

for ii = 1:ng
    
    [sdr1,sdr2] = mij2sdr(M{ii},[1 2 3]);
    sdr(ii,:)  = [sdr1 sdr2];
    crsX(:,ii) = M{ii}(indx);
    
    [~,~,xyzP(:,ii),xyzT(:,ii)] = sdr2PTaxes(sdr(ii,:));
    xyzI(:,ii) = cross(xyzP(:,ii),xyzT(:,ii));
    
    
end

jnz = find(abs(xyzI(3,:)) > 0);
xyzI(:,jnz) = -xyzI(:,jnz)./repmat(sign(xyzI(3,ii)),1,length(jnz));

save('coarseDCgrid_uniform.mat','crsX','ng','xyzP','xyzT','xyzI');
