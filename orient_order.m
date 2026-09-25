function [s,phi] = orient_order(xq,yq)
% This function defines the neamtic order field 

% paramters to define the defect core location and axis 
rotAng = 0; transX = 0.0; transY = 0.0;

% define polar grid
theta_q = atan2(yq-transY,xq-transX)+rotAng; 
rq = sqrt((xq-transX).^2+(yq-transY).^2);

K = 0.005; % defect core size
s = (1-exp(-(rq.^2)/(2*K))); % order parameters
n = -0.5; % topological charge 

% Defect axis rotation to align with experimental comparison 
if n == -0.5
    phi0 = pi/4; phi = (n*theta_q+phi0);
else 
    phi0 = -pi/4; phi = (n*theta_q+phi0);
end 

end 