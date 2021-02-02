function [wp,wsv] = vertical_P_SV_coefficients(RA,q)
% function [wp,wsv] = vertical_P_SV_coefficients(RA,q)
% 
% 2021-02-02
% This function comes from Gutenberg 1944 --- it returns coeffiecients
% to convert measured P and S amplitudes on a vertical channel to full
% P, SV amplitudes. These are typically unstable for angles between 
% 28 to 40 degrees or so (depending of Vp/Vs), near the critical angle
%
% INPUTS
%
%     RA == incidence angle at surface (0 = vertically up)
%      q == Vp/Vs ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- For incident P-wave 
a = RA;
B = asind(sind(a)/q);

% -- Equations 7b,29,28
L  = q*sqrt(cotd(2*B).*cosd(2*B)./sind(2*a));
G  = -q*(2./(L.^2+1)).*cotd(2*B);
wp = abs(G./(2*sind(B)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Now for upcoming SV wave
B = RA;
a = asind(sind(B)*q);

% -- Equations 7b, 31,30 
L = q*sqrt(cotd(2*B).*cosd(2*B)./sind(2*a));
W   = (1/q)*((2*L.^2)./(L.^2+1)).*tand(2*B);
wsv = abs(W.*(cosd(a)./cosd(2*B)));
