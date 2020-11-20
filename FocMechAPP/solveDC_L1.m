function [sdr1,sdr2,errPHI] = solveDC_L1(X,Prob,gP,gT,gI)
% function [sdr1,sdr2,errPHI] = solveDC_L1(X,Prob,gP,gT,gI)
%
% 2020-07-20
% This function finds the focal mechanism that minimizes the L1
% error angle (equivalent to the Laplace statistic "b"). This requires
% having a pre-computed grid of focal mechanisms, and the probability of
% each of them already computed. See createSDRgrid...
% 
%  INPUTS
%
%      X = initial estimated principal axes (Pressure,Tension,Interme.)
%          as columns of a 3x3 matrix (rows are x,y,z)
%   Prob = probability of a grid of NG focal mechanisms corresponding to the
%          principal stress directions gP,gT,gI
%     gP = 3 x NG matrix of pressure axes [x; y; z]
%     gT = 3 x NG matrix of tension axes [x; y; z]
%     gI = 3 x NG matrix of intermediate axes [x; y; z]
%
% OUTPUTS
%
%   sdr1,sdr2 = strike/dip/rakes of focal mechanism that produces the 
%               lowest errPHI
%      errPHI = error assigned to focal mechanism, in degrees. Equivalent
%               to Laplacian "b" statistic

% -- Largest rotation I'll try is 1 degree
dp = 2.^[4:-1:-4];
dp = [dp; -dp];
dp = dp(:);
Ndp = length(dp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- Solve initial error
ptb = abs([X(:,1)'*gP
           X(:,2)'*gT
           X(:,3)'*gI]); 

% -- It doesn't seem important to worry about angles > 90 deg       
phi = min(acosd(0.5*(sum(ptb)-1)),90);

% -- errPHI equivalent to 'b' or standard devation 
errPHI = sum(Prob.*phi);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
     
    step = 0;
    
    for ii = 1:Ndp
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RMx = [1 0 0; 0 cosd(dp(ii)) -sind(dp(ii)); 0 sind(dp(ii)) cosd(dp(ii))];
        Xt = RMx*X;
        
        ptb = abs([Xt(:,1)'*gP
                   Xt(:,2)'*gT
                   Xt(:,3)'*gI]); 

        % -- It doesn't seem important to worry about angles > 90 deg       
        phi = min(acosd(0.5*(sum(ptb)-1)),90);

        % -- errPHI equivalent to 'b' or standard devation 
        errPHIt = sum(Prob.*phi);
        
        if errPHIt < errPHI
            X = Xt;
            errPHI = errPHIt;
            step=1;
            break
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RMy = [cosd(dp(ii)) 0 sind(dp(ii)); 0 1 0; -sind(dp(ii)) 0 cosd(dp(ii))];
        Xt = RMy*X;
        
        ptb = abs([Xt(:,1)'*gP
                   Xt(:,2)'*gT
                   Xt(:,3)'*gI]); 

        % -- It doesn't seem important to worry about angles > 90 deg       
        phi = min(acosd(0.5*(sum(ptb)-1)),90);

        % -- errPHI equivalent to 'b' or standard devation 
        errPHIt = sum(Prob.*phi);
        
        if errPHIt < errPHI
            X = Xt;
            errPHI = errPHIt;
            step=1;
            break
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        RMz = [cosd(dp(ii)) -sind(dp(ii)) 0; sind(dp(ii)) cosd(dp(ii)) 0; 0 0 1];
        Xt = RMz*X;
        
        ptb = abs([Xt(:,1)'*gP
                   Xt(:,2)'*gT
                   Xt(:,3)'*gI]); 

        % -- It doesn't seem important to worry about angles > 90 deg       
        phi = min(acosd(0.5*(sum(ptb)-1)),90);

        % -- errPHI equivalent to 'b' 
        errPHIt = sum(Prob.*phi);
        
        if errPHIt < errPHI
            X = Xt;
            errPHI = errPHIt;
            step=1;
            break
        end
    end
    
    
    if ii == Ndp && step == 0
        break
    else
        dp(1:2*(fix(ii/2)-1) ) = [];
        Ndp = length(dp);
    end
    
end

% -- Now just convert P,T axes to strike/dip/rake
[sdr1,sdr2] = PTaxes2sdr(X(:,1),X(:,2));
