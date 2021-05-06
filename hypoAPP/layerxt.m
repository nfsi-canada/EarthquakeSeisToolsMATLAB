function [dx,dt,irtr] = layerxt(p,h,vtop,vbot)
% FUNCTION [DX,DT,IRTR] = LAYERXT(P,H,VTOP,VBOT)
%
% LAYERXT calculates DX, DT for a ray in a layer with a
% linear velocity gradient. This is a highly modified version
% of a routine in Chapman's WKBJ program. P is horizontal slowness,
% H is layer thickness, VTOP is velocity at top of layer, VBOT is
% velocity at bottom of layer. IRTR is return code = -1 for zero
% thickness layer, =0 for ray turning above layer, =1 for ray passed
% through layer and =2 for ray turned in layer (only 1 leg counted in 
% DX, DT).
%
% NOTE ON AUTHORSHIP:  This code was included with the resources for
% the textbook "Introduction to Seismology" by Peter M. Shearer. 

utop=1/vtop;
ubot=1/vbot;
if(p>=utop)
  dx=0.0;
  dt=0.0;
  irtr=0;
  return
elseif (h==0)
  dx=0.0;
  dt=0.0;
  irtr=-1;
  return
end

% Set key properties. 
u1=utop;
u2=ubot;
v1=vtop;
v2=vbot;
b=(v2-v1)/h;
eta1=sqrt(u1^2-p^2);

% Case of homogeneous layer (ray must 
% exit below.
if(b==0)
  dx=h*p/eta1;
  dt=h*u1^2/eta1;
  irtr=0;
end
x1=eta1/(u1*b*p);
tau1=(log((u1+eta1)/p)-eta1/u1)/b;

% Case that ray bottoms within layer (only
% count one leg).
if(p>=ubot)
  dx=x1;
  dtau=tau1;
  dt=dtau+p*dx;
  irtr=2;
  return
end

% Ray must exit through bottom of layer if
% we arrive at this point.
irtr=1;
eta2=sqrt(u2^2-p^2);
x2=eta2/(u2*b*p);
tau2=(log((u2+eta2)/p)-eta2/u2)/b;

dx=x1-x2;
dtau=tau1-tau2;
dt=dtau+p*dx;

return
end




