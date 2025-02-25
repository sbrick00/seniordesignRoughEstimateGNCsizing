% torque_magnetic.m 
%   Disturbance torque from magnetic field 
% 
%   Input: 
%     lat = vehicle latitude [rad] 
%     r = vehicle position vector magnitude [m] 
%     re = earth radius [m] 
% 
%   Output: 
%     T_mag = magnetic field torque [Nm] 
% 
function T_mag = torque_magnetic(lat,r,re) 
  %# Earth magnetic field (approx as dipole) 
    B = (1 + sin(lat)^2)^(0.5) * 0.3/((r/re)^3); %[gauss] 
    B_t = B*1e-4; %[tesla], [N/(A*m)] 
  %# Vehicle residual dipole 
    D = 2; %[A*m^2] 
  %# Mag torque 
    T_mag = B_t*D; %[Nm] 