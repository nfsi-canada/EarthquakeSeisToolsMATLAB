function Mo = M2Mo(M)
% function Mo = M2Mo(M)
%
% 2018-10-04
% This function takes in a either a 3x3 MT, a 6-element MT vector, or a 
% 5-element MT vector, and computes a scalar seismic moment. If a 6xNE or 5XNE
% matrix is provided, this can take many events at once. 
%
%     INPUTS
%     
%     either  M = [M11 M12 M13; M12 M22 M23; M13 M23 M33],
%             M = [M11 M22 M33 M12 M13 M23]', or
%             M = [M11 M22 M12 M13 M23]'
%
%     OUTPUTS
%
%            Mo = scalar seismic moment

if size(M,1) == 3
    Mo = sqrt(0.5*sum(M(:).^2));
elseif size(M,1) == 6
    Mo = sqrt(0.5*sum(M(1:3,:).^2)+sum(M(4:6,:).^2))';
elseif size(M,1) == 5
    Mo = sqrt(0.5*(sum(M(1:2,:).^2)+sum(M(1:2,:)).^2)+sum(M(3:5,:).^2))';
end


