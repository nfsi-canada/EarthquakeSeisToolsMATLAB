function [tpp,A] = RayTrace3D_Node(x0,x1,gx,gy,gz,vg,wght)
% function [tpp,A] = RayTrace3D_Node(x0,x1,gx,gy,gz,vg,wght)
%
% 2020-06-12
% This function traces rays thorugh a 3D (grid) velocity model 
% with the "bending" technique of Um and Thurber [1987,BSSA].
% It then computes a distance corresponding to each cell.
% 
% This code does not apply any earth-flattening transformation, so 
% gz and vg must be already adjusted accordingly if desired. 
%
%     INPUTS
%
%          x0 == source coordinates [x,y,z] (with Nray rows).
%          x1 == receiver coordinates [x,y,z], typically z would be 0.
%                but you can include station elevation if desired 
%                (any of x,y,z can be negative).
%    gx,gy,gz == vectors of grid points in the horizontal (x2),vertical 
%                directions
%          vg == velocity at each point of the grid defined by gx,gz
%                units should match gx,gz. Must be an "NDGRID" rather  
%                than "meshgrid", i.e. x is 1st dimension, z is 3rd.
%        wght == Nray x 1 weights to multiple equations of A by
%
%     OUTPUTS
%
%         tpp == predicted travel time (in unist matching V)
%           A == an Nray x Nnode sparse matrix assigning

% -- Create interpolant of velocity model for faster sampling later
FV = griddedInterpolant({gx,gy,gz},vg);

Nray = size(x0,1);
tpp  = Inf*ones(Nray,1);
RP   = cell(Nray,1); 

% -- Don't converge unless there are at least Nmin points in ray-path
% -- ... I could make this a function input
Nmin = 32;


% -- Pieces needed to assign ray-path to nodes.
dgx = diff(gx(1:2))/2;
dgy = diff(gy(1:2))/2;
dgz = diff(gz(1:2))/2;
Nx  = length(gx);
Ny  = length(gy);
Nz  = length(gz);
Nc  = Nx*Ny*Nz;
Mj  = repmat((1:Nx)',1,Ny,Nz) + repmat(Nx*(0:Ny-1),Nx,1,Nz)  ...
      + repmat( permute((Nx*Ny)*(0:Nz-1),[3 1 2]),Nx,Ny,1);
A = sparse(Nray,Nc);

col = @(x) x(:);

% -- Run in parallel if pool is open
% -- I think I have memory leak issues....
pool = gcp('nocreate');
if isempty(pool)
    Nworker = 0;
else
    Nworker = pool.NumWorkers;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor (ii = 1:Nray, Nworker)
%for ii = 1:Nray
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -- (1) Trace the ray (~follow Um and Thurber [1987,BSSA])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X = [x0(ii,:); (x0(ii,:)+x1(ii,:))/2; x1(ii,:)];
    V = FV(X);
  
    N     = 3; % -- number of points in ray-path (initialized at 3)
    CNVRG = 0; % -- CNVRG == 0 until path converges
    
    while ~CNVRG

        % -- Adjust points from outside in (see Um and Thurber Fig. 2)
        rngk = [2:(N-1)/2; N-1:-1:(N+3)/2];
        rngk = [rngk(:); (N+1)/2];
        iter = 1; 
        
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
                % -- step length. Conisider independent x,y,z distances
                % -- The scaling here is quite arbitrary...I'm not sure
                % -- exactly what's best
                Lx = L*max(abs(dX(1)/L),0.2)/2;
                Ly = L*max(abs(dX(2)/L),0.2)/2;
                Lz = L*max(abs(dX(3)/L),0.2)/2;
                gV = [(FV(Xm(1)+Lx/2,Xm(2),Xm(3))-FV(Xm(1)-Lx/2,Xm(2),Xm(3)))/Lx ...
                      (FV(Xm(1),Xm(2)+Ly/2,Xm(3))-FV(Xm(1),Xm(2)-Ly/2,Xm(3)))/Ly ... 
                      (FV(Xm(1),Xm(2),Xm(3)+Lz/2)-FV(Xm(1),Xm(2),Xm(3)-Lz/2))/Lz];
                   
                % -- n == (unit) pertubation direction (Um and Thurber Eq. 4)
                % -- if n0 is zero this is unstable...can't use constant 
                % -- velocity model...
                n0 = gV - dot(gV,dX)*dX/norm(dX)^2;
                n  = n0/norm(n0);

                % -- R == Pertubation distance (Um and Thurber Eq. 6)
                % -- Ensure the term added to R1 is real
                c  = (1/Vt(k+1) + 1/Vt(k-1))/2;
                R1 = (c*Vt(k)+1)/(4*c*dot(n,gV));       
                R  = -R1 + sqrt(max( R1^2 + L^2/(2*c*Vt(k)) ,0));

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
            if ttmp < tpp(ii)-1e-3
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % -- (2) Assign distance travelled through each node
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dseg = vecnorm(diff(X)')';
    Nseg = length(dseg);
    
    d_node = sparse(Nc,1);
    
     for iseg = 1:Nseg
        
        Xb = sort(X(iseg:iseg+1,:));

        % -- 
        jx = find((gx > Xb(1,1)-dgx).*(gx < Xb(2,1)+dgx));
        jy = find((gy > Xb(1,2)-dgy).*(gy < Xb(2,2)+dgy));
        jz = find((gz > Xb(1,3)-dgz).*(gz < Xb(2,3)+dgz));
        
        jn = col(Mj(jx,jy,jz));

        % -- Convert 3D grid indices to indices of gx,gy,gz vectors
        % -- This is slightly slower than referencing GX,GY,GZ grids, but
        % -- I think it's better not to have to pass large grids to workers
        jz = ceil(jn/(Nx*Ny));
        jy = ceil((jn-Nx*Ny*(jz-1))/Nx);
        jx = jn - (jy-1)*Nx - (jz-1)*Nx*Ny;

        dn2 = sum([gx(jx)-Xb(1,1),gy(jy)-Xb(1,2),gz(jz)-Xb(1,3), ... 
                   gx(jx)-Xb(2,1),gy(jy)-Xb(2,2),gz(jz)-Xb(2,3)].^2,2);
        
        wn = 1./dn2;
        
        d_node(jn) = d_node(jn) + dseg(iseg)*wn/sum(wn);        
    
    end
    A(ii,:) = wght(ii)*d_node;
    
end

