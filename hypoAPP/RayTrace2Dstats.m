function [tpp,theta,dray,thetaR,v_hyp] = RayTrace2Dstats(x0,x1,gx,gz,vg)
% function [tpp,theta,dray,thetaR,v_hyp] = RayTrace2Dstats(x0,x1,gx,gz,vg)
%
% 2020-06-12
% This function traces rays thorugh a 2D (grid) velocity model 
% with the "bending" technique of Um and Thurber [1987,BSSA]. This 
% code does not apply any earth-flattening transformation, so gz and vg
% must be already adjusted accordingly if desired. This version 
% (RayTrace2Dstats vs. simply RayTrace2D) returns take-off angles and 
% the ray-path distances, instead of just returing the ray-paths
% themselves.
%
%     INPUTS
%
%          x0 == source coordinates [x,z]
%          x1 == receiver coordinates [x,z] (typically z would be 0)
%       gx,gz == vectors of grid points in the horizontal,vertical 
%                directions
%          vg == velocity at each point of the grid defined by gx,gz
%                units should match gx,gz  
%
%     OUTPUTS
%
%         tpp == predicted travel time (in unist matching V)
%       theta == the take-off angle in degress, 0 being vertically upward
%        dray == ray-path distance in units matching the grid
%      thetaR == angle of incidence at receiver, 0 being vertically upward
%       v_hyp == P-wave velocity at the hypocenter 


% -- Create interpolant of velocity model
FV = griddedInterpolant({gx,gz},vg);

Nray   = size(x0,1);
tpp    = Inf*ones(Nray,1);
theta  = zeros(Nray,1);
dray   = zeros(Nray,1);
thetaR = zeros(Nray,1);

% -- Source velocities
v_hyp = FV(x0(:,1),x0(:,2));

Nmin = 32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:Nray
    
    X = [x0(ii,:); (x0(ii,:)+x1(ii,:))/2; x1(ii,:)];
    V = [FV(x0(ii,:)); FV((x0(ii,:)+x1(ii,:))/2); FV(x1(ii,:))];
         
    N     = 3; % -- number of points in ray-path (initialized at 3)
    CNVRG = 0; % -- CNVRG == 0 until path converges
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while ~CNVRG

        % -- Adjust points from outside in (see Um and Thurber Fig. 2)
        rngk = [2:(N-1)/2; N-1:-1:(N+3)/2];
        rngk = [rngk(:); (N+1)/2];
        iter = 1; 
        
        % -- Allow multiple iterations for a given number of points
        while 1
            
            % -- Change copies of ray-path in case travel-time increases     
            Xt = X;
            Vt = V;
            
            for k = rngk'
    
                dX = Xt(k+1,:)-Xt(k-1,:); 
                Xm = (Xt(k+1,:)+Xt(k-1,:))/2;
                L  = norm(dX)/2;
                
                % -- Compute velocity gradient at mid-point, 
                % -- but over a large distance (= 1/2 * step_length) 
                Lx = L*max(dX(1)/L,0.1);
                Lz = L*max(dX(2)/L,0.1);
                gV = [(FV(Xm(1)+Lx/2,Xm(2))-FV(Xm(1)-Lx/2,Xm(2)))/Lx ...
                      (FV(Xm(1),Xm(2)+Lz/2)-FV(Xm(1),Xm(2)-Lz/2))/Lz]; 
                
                % -- n == (unit) pertubation direction (Um and Thurber Eq. 4)
                n0 = gV - dot(gV,dX)*dX/norm(dX)^2;
                n  = n0/norm(n0);

                % -- R == Pertubation distance (Um and Thurber Eq. 6)  
                c  = (1/Vt(k+1) + 1/Vt(k-1))/2;
                R1 = (c*Vt(k)+1)/(4*c*dot(n,gV));       
                R  = -R1 + sqrt (R1^2 + L^2/(2*c*Vt(k)) );

                Xt(k,:) = Xm + R*n;
                Vt(k)   = FV(Xt(k,:));
                
            end

            % -- Trapezoidal integration for travel time (Um and Thurber Eq.2)
            di   = sqrt(sum(diff(Xt).^2,2));
            ttmp = 0.5*sum( di.*(1./Vt(1:N-1) + 1./Vt(2:N)) );

            % -- If travel-time decreased or N is small, accept step
            if ttmp < tpp(ii) || N < Nmin/2
                tpp(ii) = ttmp;
                X = Xt;
                V = Vt;
            end
            
            % -- If travel-time decreased sufficiently, iterate again with
            % -- this number of points. Otherwise either double the points,
            % -- or stop if "converged".
            if ttmp < tpp(ii)*0.9999
                iter = iter+1;
                continue
            else    
                if iter==1 && N > Nmin
                    CNVRG = 1;
                else
                    X = interp1(1:N,X,1:0.5:N);
                    V = FV(X);
                    N = 2*N-1;
                end   
                break
            end
        end
    end
 
    % -- Compute Ray-Path distance and take-off angles
    dray(ii)   = sum(di);
    theta(ii)  = 90+atan2d(X(2,2)-X(1,2),abs(X(2,1)-X(1,1)));
    thetaR(ii) = 90+atan2d(X(end,2)-X(end-1,2),abs(X(end,1)-X(end-1,1)));
    keyboard
end


% -- Note that "gradient" expects a meshgrid rather than ndgrid
% -- to we need to transpose, then transpose the result.
% -- Storing these gradients requires a lot of extra memory...maybe
% -- I need a good way to quickly get gradient from vg, FV
%[GX,GZ] = gradient(vg',diff(gx(1:2)),diff(gz(1:2)));
%FGX = griddedInterpolant({gx,gz},GX');
%FGZ = griddedInterpolant({gx,gz},GZ');

%sclGx = 0.5*diff(gx(1:2));
%sclGz = 0.5*diff(gz(1:2));
%Dx = 0.5*sclGx*[-1; 1];
%Dz = 0.5*sclGz*[-1; 1];
%Z = [0; 0];

% -- gVm == local velocity gradient at mid-point. Note I am a
%        applying (grad V)_mid in Eq. 4 and Eq. 6...I'm not sure
%        if the (grad V) in Um and Thurber Eq. 4 is supposed to be
%        computed differntly...this seems to work well though.
%gVm = [FGX(Xm), FGZ(Xm)];
%gVm = [diff(FV(Xm(1)+D,Xm(2)+Z)) diff(FV(Xm(1)+Z,Xm(2)+D))];
%gVm = [diff(FV(Xm(1)+Dx,Xm(2)+Z ))/sclGx ...
%       diff(FV(Xm(1)+Z, Xm(2)+Dz))/sclGz];