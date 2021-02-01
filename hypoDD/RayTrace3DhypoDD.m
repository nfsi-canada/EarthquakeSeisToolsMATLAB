function [tpp,g,v0] = RayTrace3DhypoDD(x0,x1,FV,minV,maxV,minZ)
% function [tpp,g,v0] = RayTrace3DhypoDD(x0,x1,FV,minV,maxV,minZ)
%
% 2020-06-12
% This function traces rays thorugh a 3D velocity model 
% with the "bending" technique of Um and Thurber [1987,BSSA]. This 
% code does not apply any earth-flattening transformation, so gz and vg
% must be already adjusted accordingly if desired. 
%
% 2021-01-06
% This function is the same as RayTrace3D (as of 2021-01-06)
% except it only returns the ray takeoff direction instead of the full
% raypath...for efficient use in hypoDD. It also returns the velocity at
% the source.
%
%     INPUTS
%
%          x0 == source coordinates [x,y,z] (with Nray rows).
%          x1 == receiver coordinates [x,y,z], z would often be 0.
%                but you can add station depth
%          FV == Velocity Model, as a gridded interpolant. I suggest 
%                making it with a command like:
%                FV = griddedInterpolant({gx,gy,gz},vg,'nearest');
%                where gx,gy,gz == vectors of grid points, and
%                vg == velocity on the grid defined by gx,gy,gz
%                units should match gx,gy,gz. Must be an "NDGRID" rather  
%                than "meshgrid", i.e. x is 1st dimension, z is 3rd.
%        minV == ensure interpolated velocity is greater than minV
%        maxV == ensure interpolated velocity is less than maxV
%        minZ == ensure no part of raypath is shallower than minZ
%                default is the station depth
%
%     OUTPUTS
%
%         tpp == predicted travel time (in units matching V)
%           g == take-off directions (Nray x 3 matrix), 
%          v0 == velocity at sources

% -- Create interpolant of velocity model for faster sampling later


Nray = size(x0,1);
tpp  = Inf*ones(Nray,1);
g    = zeros(Nray,3);
v0   = zeros(Nray,1);

% -- Default minimum velocities and minimum depth
if nargin < 6
    minZ = min(x0(:,3),x1(:,3));
    if nargin < 5
        maxV = 10;
        if nargin < 4
            minV = 3;
        end
    end
end


% -- Don't converge unless there are at least Nmin points in ray-path
% -- ... I could make this a function input
Nmin = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:Nray    

    N     = 3;
    CNVRG = 0; % -- CNVRG == 0 until path converges
    swNaN = 0;
   
    wX  = repmat(linspace(0,1,N)',1,N);
    X = flipud(wX).*repmat(x0(ii,:),N,1)+wX.*repmat(x1(ii,:),N,1);
    V = Vfun(X,FV,minV,maxV);
    
    N  = size(X,1);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while ~CNVRG

        % -- Adjust points from outside in (see Um and Thurber Fig. 2)
        rngk = [2:(N-1)/2; N-1:-1:(N+3)/2];
        rngk = [rngk(:); (N+1)/2];
                
        iter  = 1; 
        
        % -- Allow multiple iterations for a given number of points
        while 1
            
            % -- Keep a copy if pre-iteration ray-path in case T increases     
            Xt = X;
            Vt = V;
            
            for k = rngk'
    
                dX = Xt(k+1,:)-Xt(k-1,:); 
                Xm = (Xt(k+1,:)+Xt(k-1,:))/2;
                L  = norm(dX)/2;   

                % -- Compute velocity gradient at mid-point, 
                % -- but over a window that is a substantial portion of the
                % -- step length. Consider independent x,y,z distances
                % -- The scaling here is quite arbitrary...I'm not sure
                % -- exactly what's best
                LdX = 0.9*abs(dX);
                VgV =  Vfun([Xm(1)+LdX(1)/2,Xm(2),Xm(3)            
                             Xm(1),Xm(2)+LdX(2)/2,Xm(3)
                             Xm(1),Xm(2),Xm(3)+LdX(3)/2
                             Xm(1)-LdX(1)/2,Xm(2),Xm(3) 
                             Xm(1),Xm(2)-LdX(2)/2,Xm(3)
                             Xm(1),Xm(2),Xm(3)-LdX(3)/2],FV,minV,maxV);                
                gV = ((VgV(1:3)-VgV(4:6))')./LdX;
                
                % -- n == (unit) pertubation direction (Um and Thurber Eq. 4)
                n0 = gV - dot(gV,dX)*dX/norm(dX)^2;
                n  = n0/norm(n0);

                % -- R == Pertubation distance (Um and Thurber Eq. 6)
                c  = (1/Vt(k+1) + 1/Vt(k-1))/2;
                R1 = (c*Vt(k)+1)/(4*c*dot(n,gV));       
                Rn = (-R1 + sqrt(R1^2 + L^2/(2*c*Vt(k))))*n;
                
                % -- New way to deal with gradient==0 NaN issues...
                Rn(isnan(Rn)) = 0;
                
                %Rn = real(Rn);
                if ~isreal(Rn)
                   tpp(ii) = NaN;
                   swNaN   = 1;
                   break
                end
                
                Xt(k,:) = Xm + Rn;   
                Xt(k,3) = max(Xt(k,3),minZ(ii));
                Vt(k)   = Vfun(Xt(k,:),FV,minV,maxV);
                
            end

            if swNaN
                break
            end
            
            % -- Trapezoidal integration for travel time (Um and Thurber Eq.2)
            di   = sqrt(sum(diff(Xt).^2,2));
            ttmp = 0.5*sum( di.*(1./Vt(1:N-1) + 1./Vt(2:N)) );

            % -- If travel-time decreased or N is small, accept step
            % -- Sharp velocity gradients may result in travel-time
            % -- increases, its better to accept anyway for small N.
            if ttmp < tpp(ii) || N < Nmin/2
                tpp(ii) = ttmp;
                X = Xt;
                V = Vt;
            end
            
            % -- If travel-time decreased sufficiently, iterate again with
            % -- this number of points. Otherwise either double the points,
            % -- or stop if "converged".
            if ttmp < tpp(ii)*0.999
                iter = iter+1;
                continue
            else    
                if iter==1 && N > Nmin
                    CNVRG = 1;
                else
                    X = interp1(1:N,X,1:0.5:N);
                    V = Vfun(X,FV,minV,maxV);
                    N = 2*N-1;
                end   
                break
            end
        end
        
        if swNaN
            break
        end
        
    end
    
    % -- Save take-off angle
    %RP{ii} = X;
    dX = diff(X(1:2,:));
    g(ii,:) = dX/norm(dX);
    v0(ii)  = V(1);
end





