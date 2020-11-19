function [sdrs1,sdrs2,vpvs] = mij2sdrs(M,cnvtn)
% function [sdrs1,sdrs2,vpvs] = mijsdrs(M,cnvtn)
%
% This function calculates a strike, dip, and rake from a moment tensor.
% It also computes the slope, which is the angle of the
% dislocation vector from the fault plane. It also takes in a 3-integer 
% vector describing the order and polarity of the coordinate system used.
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
%       sdrs1,2 = the complementary solutions [strk,dip,rke,slp]
%
%       str = fault strike using right-hand-rule (0-360 clockwise from N)
%       dip = fault dip using right ahnd rule (0-90 from horizontal)
%       rke = fault rake, describes sense of motion, describes motion of upper
%             block relative to strike. Clockwise, in degrees. 
%                 0       = left-lateral
%                 +/- 180 = right-lateral
%                 +90     = reverse
%                 -90     = normal
%       slp = slope, in degrees, angle from fault plane of dislocation vector
%                 0 = shear, +90 = extensional, -90 = compressional
%
%      vpvs = the Vp/Vs ratio implied by this moment-tensor





% THIS ISN'T WORKING FOR  |RAKES| > 90! ...or maybe sdrs2mij isn't? 
% Not working for rakes < -90, for rakes == 180

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

[V,D] = eig(M);

[~,jj] = sort(diag(D),'descend');
V = V(:,jj);
D = D(jj,jj);


M1 = D(1,1);
M2 = D(2,2);
M3 = D(3,3);

% -- Check if consistent with tensile source model
% -- Not sure if 
%c = abs((M1+M3-2*M2)/(M1-M3))*(M1+M2+M3)/(M1+M3-2*M2);

%if c < 0
%    keyboard
%end


% -- Get Vp/Vs ratio
vpvs = sqrt(1 + (M1+M3)/(M1+M3-2*M2));

% -- Get the slope
slp = asind((M1+M3-2*M2)/(M1-M3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Get first solution

% -- From Vavrycuk 2011 (n == fault normal, v = dislocation direction)
n = sqrt((M1-M2)/(M1-M3))*V(:,1) + sqrt((M3-M2)/(M3-M1))*V(:,3);
v = sqrt((M1-M2)/(M1-M3))*V(:,1) - sqrt((M3-M2)/(M3-M1))*V(:,3);

% -- First convert to E,N,U if necessary
%[~,jj] = ismember([1,2,3],abs(cnvtn));
%sgn    = sign(cnvtn(:));
%n      = sgn(jj).*n(jj);
%v      = sgn(jj).*v(jj);

% -- Ensure normal is pointing upward 
if n(3) < 0
    n = -n;
    v = -v;  % -- Not sure this is a good approach...
end

% -- Dip can now be taken without complication
dip = acosd(n(3));

% -- Horizontal angle of normal (dip direction) from -180 to 180 (W to W) counter-clockwise
phi = atan2d(n(2),n(1));
phi = 90-phi; % -- dip direction in map angle, from -90 to 270 (W to W) clockwise
str = phi-90; % -- Strike from dip direction (-180 to 180) (S to S) clockwise
if str < 0
    str = 360+str;
end

% -- Rake
rke = asind((v(3)-cosd(dip)*sind(slp))/(sind(dip)*cosd(slp)));
cosR = v(1)*cosd(90-str)+v(2)*sind(90-str);
rke(cosR<0) = sign(rke)*180-rke;
%{
% -- Correct rake if necessary (convert to |rke| > 90)
rphi = atan2d(v(2),v(1));
rphi = 90-rphi;

rphi(rphi < 0)      = 360+rphi;
rrphi               = rphi-str;
rrphi(rrphi >  180) = -360+rrphi; 
rrphi(rrphi < -180) =  360+rrphi; 

%rphi(rphi < 0)      = 360+rphi;
%rrphi               = rphi-phi;
%rrphi(rrphi >  180) = -360+rrphi; 
%rrphi(rrphi < -180) =  360+rrphi; 

if rke > 0 && rrphi < -90
    rke = 180-rke;
elseif rke < 0 && rrphi > 90
    rke = -180-rke;
end
%}

sdrs1 = [str,dip,rke,slp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Get second solution (swap n and v)

% -- From Vavrycuk 2011 (n == fault normal, v = dislocation direction)
v = sqrt((M1-M2)/(M1-M3))*V(:,1) + sqrt((M3-M2)/(M3-M1))*V(:,3);
n = sqrt((M1-M2)/(M1-M3))*V(:,1) - sqrt((M3-M2)/(M3-M1))*V(:,3);

% -- First convert to E,N,U if necessary
%[~,jj] = ismember([1,2,3],abs(cnvtn));
%sgn    = sign(cnvtn(:));
%n      = sgn(jj).*n(jj);
%v      = sgn(jj).*v(jj);

% -- Ensure normal is pointing upward 
if n(3) < 0
    n = -n;
    v = -v; 
end

% -- Dip can now be taken without complication
dip = acosd(n(3));

% -- Horizontal angle of normal (dip direction) from -180 to 180 (W to W) counter-clockwise
phi = atan2d(n(2),n(1));
phi = 90-phi; % -- dip direction in map angle, from -90 to 270 (W to W) clockwise
str = phi-90; % -- Strike from dip direction (-180 to 180) (S to S) clockwise
str(str < 0) = 360+str;


% -- Rake, This is ambiguious--if |rke| > 90 it will be wrong (sin symmetric around 90)
rke = asind((v(3)-cosd(dip)*sind(slp))/(sind(dip)*cosd(slp)));
cosR = v(1)*cosd(90-str)+v(2)*sind(90-str);
rke(cosR<0) = sign(rke)*180-rke;
%{
% -- Correct rake if necessary (convert to |rke| > 90)
rphi = atan2d(v(2),v(1));
rphi = 90-rphi;

rphi(rphi < 0)      = 360+rphi;
rrphi               = rphi-str;
rrphi(rrphi >  180) = -360+rrphi; 
rrphi(rrphi < -180) =  360+rrphi; 

%rphi(rphi < 0)      = 360+rphi;
%rrphi               = rphi-phi;
%rrphi(rrphi >  180) = -360+rrphi; 
%rrphi(rrphi < -180) =  360+rrphi; 

if rke > 0 && rrphi < -90
    rke = 180-rke;
elseif rke < 0 && rrphi > 90
    rke = -180-rke;
end
%}

sdrs2 = [str,dip,rke,slp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% -- For consitency, always give solution with lower strike as sdrs1
if sdrs1(1) > sdrs2(1)
    tmp   = sdrs1;
    sdrs1 = sdrs2;
    sdrs2 = tmp;
end

