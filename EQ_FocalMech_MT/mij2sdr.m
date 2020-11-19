function [sdr1,sdr2] = mij2sdr(M,cnvtn)
% function [sdrs1,sdrs2,vpvs] = mijsdr(M,cnvtn)
%
% This function calculates a best fit DC from a moment tensor.
% It also takes in a 3-integer vector describing the order and polarity 
% of the coordinate system used.
%
%    INPUTS
%
%         M = the moment tensor
%     cnvtn = coordinate order and polarity, relative to E,N,U
%             E,N,U == [1,2,3]
%             N,E,D == [2,1,-3]
%
%    OUTPUTS
%
%       sdrs1,2 = the complementary solutions [strk,dip,rke]
%
%       str = fault strike using right-hand-rule (0-360 clockwise from N)
%       dip = fault dip using right ahnd rule (0-90 from horizontal)
%       rke = fault rake, describes sense of motion, describes motion of upper
%             block relative to strike. Clockwise, in degrees. 
%                 0       = left-lateral
%                 +/- 180 = right-lateral
%                 +90     = reverse
%                 -90     = normal

% -- Convert to ENU if necessary
if nargin > 2
    [~,jj]  = ismember([1,2,3],abs(cnvtn));
    pj      = sign(cnvtn(jj));
    rm      = zeros(3);
    for ii = 1:3
        rm(ii,jj(ii)) = pj(ii);
    end
    M = rm*M*rm';
end

% -- Remove volumetric part, take eigenvalues
M     = M - eye(3)*trace(M)/3;
[V,D] = eig(M);


[~,jj] = sort(diag(D),'descend');
V = V(:,jj);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- 1st Plane

% -- From Vavrycuk 2011 (n == fault normal, v = dislocation direction)
% -- Modified to assume M2 = 0, M1 = -M3
n = (V(:,1) + V(:,3))/sqrt(2);
v = (V(:,1) - V(:,3))/sqrt(2);

% -- Ensure normal is pointing upward 
if n(3) < 0
    n = -n;
    v = -v;  
end

% -- Dip can now be taken without complication
dip = acosd(n(3));

% -- Horizontal angle of normal (dip direction) from -180 to 180 (W to W) counter-clockwise
str = -atan2d(n(2),n(1));
str(str<0) = str+360;

rke  = asind(v(3)/sind(dip));
cosR = v(1)*cosd(90-str)+v(2)*sind(90-str);
rke(cosR<0) = sign(rke)*180-rke;
%rke(cosR<0)   = 180-rke;
%rke(rke<-180) = rke+360;
%rke(rke>180)  = rke-360;

sdr1 = [str,dip,rke];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Get second solution (swap n and v)


% -- From Vavrycuk 2011 (n == fault normal, v = dislocation direction)
% -- Modified to assume M2 = 0, M1 = -M3
n = (V(:,1) - V(:,3))/sqrt(2);
v = (V(:,1) + V(:,3))/sqrt(2);


% -- Ensure normal is pointing upward 
if n(3) < 0
    n = -n;
    v = -v;  
end

% -- Dip can now be taken without complication
dip = acosd(n(3));

% -- Horizontal angle of normal (dip direction) from -180 to 180 (W to W) counter-clockwise
str = -atan2d(n(2),n(1));
str(str<0) = str+360;

% -- Rake
rke  = asind(v(3)/sind(dip));
cosR = v(1)*cosd(90-str)+v(2)*sind(90-str);
rke(cosR<0) = sign(rke)*180-rke;
%rke(cosR<0)   = 180-rke;
%rke(rke<-180) = rke+360;
%rke(rke>180)  = rke-360;

sdr2 = [str,dip,rke];


