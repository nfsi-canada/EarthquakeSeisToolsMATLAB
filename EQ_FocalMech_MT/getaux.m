%% This code is taken from Shearer, but translated from Fortran into Matlab
%% GETAUX returns auxiliary fault plane, given strike,dip,rake of main fault plane.  
%% Inputs:  strike1, dip1, rake1 (degrees, primary fault plane)
%% Returns: strike2, dip2, rake2 (degrees, auxiliary fault plane)
%% Strike, using right-hand-rule
%% Rake: 0=left lateral, 90=reverse, +-180=right lateral, - 90=normal
%% If only one (three-element) 

function [s2,d2,r2] = getaux(s1,d1,r1)

    if nargin == 1
        r1 = s1(:,3);
        d1 = s1(:,2);
        s1 = s1(:,1);
    end

    d2  = acosd(sind(r1).*sind(d1));
   
    sr2 = cosd(d1)./sind(d2);
    cr2 = -sind(d1).*cosd(r1)./sind(d2);
    r2  = atan2d(sr2,cr2);
   
    s12 = cosd(r1)./sind(d2);
    c12 = -1./(tand(d1).*tand(d2));
    s2  = s1-atan2d(s12,c12);

    for ii = 1:length(d2)
        if d2(ii) > 90 
            s2(ii) = s2(ii)+180;
            d2(ii) = 180-d2(ii);
            r2(ii) = 360-r2(ii);
        end

        if s2(ii) > 360
            s2(ii) = s2(ii)-360;
        end

        if r2(ii) > 180
            r2(ii) = r2(ii)-360;
        end
    end

    if nargin == 1
        s2 = [s2, d2, r2];
    end
end
%{

ORIGINAL FORTRAN CODE
      subroutine GETAUX(strike1,dip1,rake1,strike2,dip2,rake2)
      degrad=180./3.1415927
      s1=strike1/degrad
      d1=dip1/degrad
      r1=rake1/degrad

      d2=acos(sin(r1)*sin(d1))

      sr2=cos(d1)/sin(d2)
      cr2=-sin(d1)*cos(r1)/sin(d2)
      r2=atan2(sr2,cr2)

      s12=cos(r1)/sin(d2)
      c12=-1./(tan(d1)*tan(d2))
      s2=s1-atan2(s12,c12)

      strike2=s2*degrad
      dip2=d2*degrad
      rake2=r2*degrad

      if (dip2.gt.90.) then
         strike2=strike2+180.
         dip2=180.-dip2
         rake2=360.-rake2
      end if
      if (strike2.gt.360.) strike2=strike2-360.

      return
      end 
%}

