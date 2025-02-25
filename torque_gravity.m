% torque_gravity.m 
%   Disturbance torque from gravity gradient 
% 
%   Input: 
%     r = vehicle radius [m] 
%     mu = planet gravity constant [m^3/s^2] 
%     I = vehicle moment of inertia [kg*m^2] 
% 
%   Output: 
%     T_grav = gravity gradient torque [Nm] 
% 
function T_grav = torque_gravity(r, mu, I) 
  %# Max moment 
    Imax = max(diag(I)); %[kg*m^2] 
  %# Min moment 
    Imin = min(diag(I)); %[kg*m^2] 
  %# Angle deviation from vertical 
    theta = 45*pi/180; %[rad], worst case angle chosen 
  %# Calc gravity gradient torque 
    T_grav = 3*mu*sin(2*theta)*(Imax - Imin)/(2*r^3); %[N*m]