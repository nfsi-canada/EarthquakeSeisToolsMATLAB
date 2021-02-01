function [tpp,theta,dray,thetaR,v_hyp] = RayTrace1D(hdepth,xep,model)
% function [tpp,theta,dray,thetaR,v_hyp] = RayTrace1D(hdepth,xep,model)
%
% 2020-06-10
% This function traces rays thorugh a layered 1D (gradient) velocity
% model with the "bending" technique of Um and Thurber [1987,BSSA]. 
% This is a more robust method than used in my "RayTrace" code, but much 
% less efficient if you need multiple rays computed. This applies an
% earth-flattening transformation.
%
% I should alter this to take in x0,x1 like RayTrace2D
%
%   INPUTS
%
%    hdepth == hypocentral depth in km
%       xep == epicentral distance in km
%     model == name of velocity model e.g. 'GSC','ALASKA','ma2011'
%                OR a custom gradient model [Depth (km), Velocity (km s-1)]
%   OUTPUTS
%
%       tpp == predicted P travel time in s
%        RP == ray-path: a list of position [X (from origin),Z]

if ischar(model)
    switch model
    case 'ma2011'
    DV = [  0.0  6.25
            8.0  6.50
           17.0  6.60
           24.0  6.70
           30.0  7.10
           36.0  7.20
           37.0  8.00
           60.0  8.05
          185.0  8.10
          220.0  8.55
          600.0 10.15];           
    case 'GSC'
        DV = [0.0  5.00
          1.0  5.50
          6.0  6.35
         30.0  6.90
         45.0  7.42
         50.0  7.77
         60.0  7.78
        100.0  8.00
        600.0 12.00];
    case 'ALASKA'
        DV = [0.00  5.8
          24.3  6.0
          24.4  7.3
          40.1  7.5
          40.2  7.8
          75.9  8.0
          76.0  8.27
         300.0  8.31
         301.0  10.3
         524.0  10.5       
         525.00 12.6
        1000.00 13.0];
    case 'customLSL'
    DV = [  0.0  6.00
            2.0  6.10
            8.0  6.50
           17.0  6.70
           24.0  6.80
           30.0  7.10
           36.0  7.30
           40.0  8.00
           60.0  8.05
          185.0  8.10
          220.0  8.55
          600.0 10.15];   
    end 
else
    DV = model;
end
depth = DV(:,1);
vp    = DV(:,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Flatten earth model
a   = 6371; % km
rf  = a-depth;
zf  = -a*log(rf/a);
vf  = (a./rf).*vp;
hzf = -a*log((a-hdepth)/a);

% -- Create functions for velocity, gradient
FV = griddedInterpolant(zf,vf);

sclG = 0.1*diff(zf(1:2));
D = 0.5*sclG*[-1; 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nray   = length(xep(:));
tpp    = Inf*ones(size(xep));

RP = cell(Nray,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:Nray
    
    X    = [xep(ii) hzf(ii); xep(ii)/2 hzf(ii)/2; 0 0];
    V    = [FV(hzf(ii)); FV(hzf(ii)/2); FV(0)];
    
    N    = 3;

    CNVRG = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while ~CNVRG

        % -- Adjust points from outside in (see Um and Thurber Fig. 2)
        rngk = [2:(N-1)/2; N-1:-1:(N+3)/2];
        rngk = [rngk(:); (N+1)/2];
        iter = 1;
        
        while 1
            
            X0 = X;
            V0 = V;
            
            for k = rngk'

                dX = X(k+1,:)-X(k-1,:); 
                Xm = (X(k+1,:)+X(k-1,:))/2;
                % -- gVm == local velocity gradient at mid-point. Note I am a
                %        applying (grad V)_mid in Eq. 4 and Eq. 6...I'm not sure
                %        if the (grad V) in Um and Thurber Eq. 4 is supposed to be
                %        computed differntly...this seems to work well though.
                %        Horizontal velocity gradient is zero in 1D model. 
                %gVm = [0,FGZ(Xm(2))];
                gVm = [0,diff(FV(Xm(2)+D))/sclG];
                
                % -- n == (unit) pertubation direction (Um and Thurber Eq. 4)
                n0 = gVm - dot(gVm,dX)*dX/norm(dX)^2;
                %n0 = gV - dot(gV,dX)*dX/norm(dX)^2;
                n  = n0/norm(n0);

                % -- THERE COULD BE A BUG!!!! I'm finding that R~0.8
                %    is better, but Um and Thurber suggest INCREASING R
                %    ...it often "converges" after 1-2 steps

                % -- R == Pertubation distance (Um and Thurber Eq. 6)
                L  = norm(X(k+1,:)-X(k-1,:))/2;
                c  = (1/V(k+1) + 1/V(k-1))/2;
                
                R1 = (c*V(k)+1)/(4*c*dot(n,gVm));       
                R  = -R1 + sqrt (R1^2 + L^2/(2*c*V(k)) );

                X(k,:) = Xm + R*n;
                V(k)   = FV(X(k,2));

            end

            % -- Trapezoidal integration for travel time (Um and Thurber Eq.2)
            di   = sqrt(sum(diff(X).^2,2));
            ttmp = 0.5*sum( di.*(1./V(1:N-1) + 1./V(2:N)) ) ;
            
            % -- If travel-time decreased sufficiently, continue
            % -- If travel-time decreased a tiny bit, accept step and stop
            % -- If travel-time increased, reset step and stop
            if ttmp < tpp(ii)*0.9999
                tpp(ii) = ttmp;
                iter = iter+1;
                continue    
            elseif ttmp > tpp(ii)
                X = X0;
                V = V0;
                if iter==1 && N > 17
                    CNVRG = 1;
                end
                break
            else
                tpp(ii) = ttmp;
                if iter==1 && N > 17
                    CNVRG = 1;
                end
                break
            end
             
        end 
        
        if ~CNVRG
            X = interp1(1:N,X,1:0.5:N);
            V = FV(X(:,2));
            N = 2*N-1;
        end
        
    end
    
    % -- Save ray-paths
    RP{ii} = X;

end
