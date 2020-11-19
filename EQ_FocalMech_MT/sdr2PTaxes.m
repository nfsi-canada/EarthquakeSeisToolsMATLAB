function [P,T,Px,Tx] = sdr2PTaxes(sdr)
% function [P,T] = sdr2PTaxes(sdr)
%
% 2020-06-03
% This code computes (DC focal mechanism) pressure and tension axes 
% from a strike/dip/rake. It always gives the downward pointing vector.
%
%    INPUT
%
%      sdr = [strike,dip,rake] OR [s1,d1,r1,s2,d2,r2]
%
%   OUTPUTS
%
%      P = [trend,plunge] of pressure axis 
%      T = [trend,plunge] of tension axis
%     Px = [x,y,z] of pressure axis
%     Tx = [x,y,z] of tension axis

Ne = size(sdr,1);

if size(sdr,2)==3
    sdr0 = sdr;
    sdr = [sdr sdr];
    for ie = 1:Ne
        M = sdrs2mij([sdr0(ie,:) 0]);
        [sdr1,sdr2] = mij2sdr(M,[1 2 3]);
        sdr(ie,:) = [sdr1 sdr2];
    end
end

no1 = [cosd(-sdr(:,1)).*sind(sdr(:,2)), sind(-sdr(:,1)).*sind(sdr(:,2)), cosd(sdr(:,2))];
no2 = [cosd(-sdr(:,4)).*sind(sdr(:,5)), sind(-sdr(:,4)).*sind(sdr(:,5)), cosd(sdr(:,5))];

% -- NEED TO CORRECTLY IDENTIFY P AND T AXES!!!
% -- First assume normal faulting
Px = -(no1+no2)/sqrt(2);
Tx = (no1-no2)/sqrt(2);
fT = -sign(Tx(:,3));
fT(fT==0) = 1;
Tx = Tx.*repmat(fT,1,3);

% -- Then find indices of reverse faulting, and swap P-T
jj = find(sdr(:,3) > 0);
tmpP = Px(jj,:);
Px(jj,:) = Tx(jj,:);
Tx(jj,:) = tmpP;

P = [90-atan2d(Px(:,2),Px(:,1)),asind(-Px(:,3))];
T = [90-atan2d(Tx(:,2),Tx(:,1)),asind(-Tx(:,3))];

P(P(:,1)<0,1) = P(P(:,1)<0,1) + 360;
T(T(:,1)<0,1) = T(T(:,1)<0,1) + 360;

   