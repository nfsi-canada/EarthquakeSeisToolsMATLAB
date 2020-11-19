function x = AngelierStressSolver(n,v,x,minbA)
% function x = AngelierStressSolver(n,v,x,minbA)
%
% 2020-08-07
% Trying to maximize "S". Based on the misfit decribed by Angelier [2002],
% but using a iterative grid-search approach
%
%
n1 = n(:,1);
n2 = n(:,2);
n3 = n(:,3);
v1 = v(:,1);
v2 = v(:,2);
v3 = v(:,3);

A = sum([n1.*v1-n3.*v3, n1.*v2+n2.*v1, n1.*v3+n3.*v1, n2.*v2-n3.*v3, n2.*v3+n3.*v2]);

X = x([1 2 3; 2 4 5; 3 5 1]);
X(3,3) = -x(1)-x(4);

[V,D] = eig(X);
D = flip(diag(D));
D = D/norm(D);
V = fliplr(V);

% -- Compute initial sum
X = V*diag(D)*V';
x = X([1 2 3 5 6]');
S = A*x;

ngR = 11;
ngA = 11;
bR = 0.5;
bA = 90;

xg = zeros(5,ngR*ngA^3);


Nstep =0;
while 1
    
    jj = 1;
      
    % -- Create 4D grid of rotations + change in R
    ra = linspace(-bA,bA,ngA);
    
    % -- How to change R? Keep sum of trace constant...sum of squares constant
    Dg = repmat(D,1,ngR);
    Dg(2,:) = linspace(max(D(3),D(2)-bR),min(D(1),D(2)+bR),ngR);
    Dg = Dg-repmat(mean(Dg),3,1);
    Dg = Dg./repmat(vecnorm(Dg),3,1);

    for iA1 = 1:ngA  
        RM1 = [1             0            0
               0 cosd(ra(iA1)) -sind(ra(iA1))
               0 sind(ra(iA1))  cosd(ra(iA1))];
        Vt1 = RM1*V;

        for iA2 = 1:ngA
            RM2 = [ cosd(ra(iA2)) 0 sind(ra(iA2))
                          0       1      0
                   -sind(ra(iA2)) 0 cosd(ra(iA2))];
            Vt2 = RM2*Vt1;

            for iA3 = 1:ngA    
                RM3 = [cosd(ra(iA3)) -sind(ra(iA3)) 0
                       sind(ra(iA3))  cosd(ra(iA3)) 0
                                0              0    1];
                Vt3 = RM3*Vt2;

                for iR = 1:ngR
                    Xt = Vt3*diag(Dg(:,iR))*Vt3';
                    xg(:,jj) = Xt([1 2 3 5 6]);  
                    jj = jj+1;
                end
            end
        end
    end
    
    Sg = A*xg;
     
    [maxSg,jmax] = max(Sg);
   
    if maxSg > S
        S = maxSg;
        x = xg(:,jmax);
        
        X = [x(1) x(2) x(3); x(2) x(4) x(5); x(3) x(5) -x(1)-x(4)];
        [V,D] = eig(X);
        
        D = flip(diag(D));
        D = D/norm(D);
        V = fliplr(V);
        
        Nstep = Nstep+1;
    elseif bA > minbA
        bA = bA/2;
        bR = 3*bR/4;
    else
        break
    end
end
Nstep             
             
             


%{

H = A'*A;


M0 = sqrt(sum(x.^2)+x(1)*x(4));

S = A*x;
b = S-1;

max_try = 8;
CNVRG = 0;
iter = 0;
L = 0.02;
while ~CNVRG

    %J  = -(H*x-A'*b);
    %dx = -H\J;
    %f1 = x'*H*dx - dx'*A'*b;
    %f2 = dx'*H*dx;
    %dx = (f1/f2)*dx;
    dx = -L*A/sqrt(sum(A.^2)+A(1)+A(4));

    CNVRG = 1;
    for ii = 1:max_try
        
        xt = x + dx*2^(1-ii);
        Mt = sqrt(sum(xt.^2)+xt(1)*xt(4));
        xt = xt*M0/Mt;
       
  
        St = A*xt;

        if St < S
            x = xt;
            S = St;
            b = S-1;
            iter = iter + 1;
            CNVRG = 0;
            break
        end
    end
end
keyboard
%M0 = S11^2 + S12^2 + S13^2 + S22^2 + S23^2 + S11*S22;
%}