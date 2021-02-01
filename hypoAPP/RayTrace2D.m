function [tpp,RP] = RayTrace2D(x0,x1,gx,gz,vg)
% function [tpp,RP] = RayTrace2D(x0,x1,gx,gz,vg)
%
% 2020-06-12
% This function traces rays thorugh a 2D (grid) velocity model 
% with the "bending" technique of Um and Thurber [1987,BSSA]. This 
% code does not apply any earth-flattening transformation, so gz and vg
% must be already adjusted accordingly if desired. 
%
%
%     INPUTS
%
%          x0 == source coordinates [x,z] (with Nray rows).
%          x1 == receiver coordinates [x,z] (typically z would be 0)
%       gx,gz == vectors of grid points in the horizontal,vertical 
%                directions
%          vg == velocity at each point of the grid defined by gx,gz
%                units should match gx,gz. Must be an "NDGRID" rather  
%                than "meshgrid", i.e. x is 1st dimension, z is 2rd.  
%
%     OUTPUTS
%
%         tpp == predicted travel time (in unist matching V)
%          RP == an Nray x 1 cell array, with each cell consiting of
%                [x,z] points of that ray-path.

% -- Create interpolant of velocity model for faster sampling later
FV = griddedInterpolant({gx,gz},vg);

Nray = size(x0,1);
tpp  = Inf*ones(Nray,1);
RP   = cell(Nray,1); 

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
                Lx = L*max(abs(dX(1)/L),0.2)/2;
                Lz = L*max(abs(dX(2)/L),0.2)/2;
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
    
    % -- Save-ray-path
    RP{ii} = X;
    
end

