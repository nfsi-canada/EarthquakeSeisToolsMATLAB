function sdr = nv2sdr(n,v)
% function sdr = nv2sdr(n,v)
%
% 2020-08-06
% Converts fault normal and rake vectors into a strike/dip/rake.
% Can take in an arbirtary number (N) of vector sets.
%
%  INPUTS
%
%     n = N x 3 matrix of fault normals [x,y,z]
%     v = N x 3 matrix of slip vectors  [x,y,z]
%         THIS CODE DOES NOT CHECK TO ENSURE THAT THE SLIP VECTORS
%         ARE IN FACT ON THE PLANE...
%
%  OUTPUTS
%
%     sdr = N x 3 matrix of [strike,dip,rake] angles, in degrees


% -- Ensure normal is pointing upward 
jn = find(n(:,3)<0);
n(jn,:) = -n(jn,:);
v(jn,:) = -v(jn,:);

% -- Dip can now be taken without complication
dip = acosd(n(:,3));

% -- Horizontal angle of normal (dip direction) 
% -- from -180 to 180 (W to W) counter-clockwise
str = -atan2d(n(:,2),n(:,1));
str(str<0) = str(str<0)+360;

rke  = asind(v(:,3)./sind(dip));
cosR = v(:,1).*cosd(90-str)+v(:,2).*sind(90-str);
rke(cosR<0) = sign(rke(cosR<0))*180-rke(cosR<0);

sdr = [str,dip,rke];
