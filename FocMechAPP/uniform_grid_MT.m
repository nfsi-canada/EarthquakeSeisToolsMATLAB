function M = uniform_grid_MT(N)
% function M = uniform_grid_MT(N)
% 
% 2020-07-21
% This function ~uniformly samples a sphere by randomly rotating a single
% initial MT (with a uniform distribution of angles around 3 axes).
%
%  INPUTS
%
%     N == number of focal mechanisms
%
%   OUTPUT
%
%    M == a cell array of moment tensors
%
% 



%R = [1-2.*qy.^2-2*qz.^2  2*qx.*qy-2*qz.*qw   2*qx.*qz+2*qy.*qw
%     2*qx.*qy+2*qz.*qw   1-2*qx.^2-2*qz.^2   2*qy.*qz-2*qx.*qw
%     2*qx.*qz-2*qy.*qw   2*qy.*qz+2*qx.*qw   1-2*qx.^2-2*qy.^2];


% -- Make sure N is even
if mod(N,2)
    N = N+1;
end

M0 = [1 0 0; 0 -1 0; 0 0 0];
M = cell(N,1);

% -- Generate (uniform) random quaternion (w,x,y,z)
q = RandomQuaternion(N/2);

for ii = 1:N/2 
        
    % -- Quaternion to rotation matrix. 
    R = Quat2RotMtrx(q(ii,:));
    
    M{2*ii-1} = R'*M0*R;
    M{2*ii}   = -M{2*ii-1};
end


% -- Algorithm from Kuffner 2004,IEEE Conf. on Robotics and Automation
    %{
    s = rand;
    s1 = sqrt(1-s);
    s2 = sqrt(s);
    t1 = 2*pi*rand;
    t2 = 2*pi*rand;
    w  = cos(t2)*s2;
    x  = sin(t1)*s1;
    y  = cos(t1)*s1;
    z  = sin(t2)*s2;
    %}
    

% -- Just use random, uniform distributions of rotation angles??
% -- Can I design an actually uniform distribution?
%a = 2*pi*rand(N,1);
%b = acos(sign(randn(N,1)).*rand(N,1));
%c = 2*pi*rand(N,1);

%R1 = [1 0 0; 0 cos(a(ii)) sin(a(ii)); 0 -sin(a(ii)) cos(a(ii))];
%R2 = [cos(b(ii)) 0 sin(b(ii)); 0 1 0; -sin(b(ii)) 0 cos(b(ii))];
%R3 = [cos(c(ii)) -sin(c(ii)) 0; sin(c(ii)) cos(c(ii)) 0; 0 0 1];

%M{ii} = R1'*(R2'*(R3'*M0*R3)*R2)*R1;


%{
% -- Make sure N is even
if mod(N,2)
    N = N+1;
end

GoldenRatio = (1 + sqrt(5))/2;

ii = [1:N]'-0.5;

tht = acos(1-2*ii/N);
phi = 2*pi*ii/GoldenRatio;

M0 = [1 0 0; 0 -1 0; 0 0 0];
M = cell(N,1);
for ii = 1:N
    R1 = [cos(phi(ii)) -sin(phi(ii)) 0; sin(phi(ii)) cos(phi(ii)) 0; 0 0 1];
    R2 = [cos(tht(ii)) 0 sin(tht(ii)); 0 1 0; -sin(tht(ii)) 0 cos(tht(ii))];
    M{ii} = R2'*(R1'*M0*R1)*R2;
end
%}
