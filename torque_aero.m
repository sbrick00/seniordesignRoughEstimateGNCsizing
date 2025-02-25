% torque_aero.m 
% Bill Nadir 
% 16.851 Satellite Engineering 
% 10/11/2003 
% Module for calculating external spacecraft torque caused by Aerodynamic forces 
function Ta = torque_aero(size,V,h,CG) 
% Here the force on the S/C, F, is calculated 
% INPUTS 
% size  = edge lengths of the S/C: vector (x,y,z) side lengths (meters) 
% V     = S/C velocity (m/s) 
% h     = S/C altitude (m) 
% CG    = location of the center of gravity (x,y,z) for the S/C (assumed offset from the 
%         geometric center of (0,0,0)) (m) 
% OUTPUT 
% Ta    = Torque on S/C due to aerodynamic drag (Nm) 
% rho is the atmospheric density at the location of the S/C 
% An atmospheric model for the upper atmosphere (h>25000m) is used to 
% approximate the density of the upper atmosphere 
% T is the atmospheric temperature, p is atmospheric pressure 
T   = -131.21 + .00299*h; % in deg C 
p   = 2.488*(((T + 273.1)/216.6)^-11.388); % pressure in KPa 
rho = p / (.2869*(T + 273.1)); % in kg/m^3 
% C_D is the drag coefficient of the cube-shaped S/C (assumed = 2.2) 
C_D = 2.2; 
% Cpa   = location of the center of aerodynamic pressure (x,y,z) 
%         (assumed at the center of the face of one side of the cube which is 
%         facing directly into the atmosphere = max drag) 
% Here the surface areas of the sides of the S/C are determined 
% This is used to assume the worst-case drag on the vehicle 
% [x*z y*z x*y] => find max 
area = [size(1)*size(3) size(2)*size(3) size(1)*size(2)]; 
max_area = max(area); 
F = 0.5*rho*C_D*(max_area^2)*(V^2); 
% here the external aerodynamic torque on the S/C is calculated 
Ta = F*max(abs(CG));