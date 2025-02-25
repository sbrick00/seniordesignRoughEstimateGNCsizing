% prop_system.m 
% Bill Nadir 
% 16.851 Satellite Engineering 
% 10/11/2003 
% Module for calculating spacecraft propulsion system mass for required 
% momentum dumping 
function [F, p_mass] = prop_system(size,CG,lifetime,H,sat_rate) 
% INPUTS 
% size      = edge lengths of the S/C: vector (x,y,z) side lengths (meters) 
% CG        = (x,y,z) coordinates of the location of the CG (offset from the 
%             geoometric center of the S/C) 
% lifetime  = required lifetime of the spacecraft [yrs] 
% H         = maximum stored momentum in any one momentum wheel (saturation 
%             point of a momentum wheel) [N*m*s] 
% sat_rate  = The rate of saturation of a momentum wheel (used to determine 
%             how often momentum needs to be dumped) [days/saturation] 
% OUTPUTS 
% p_mass    =  total mass of the propulsion subsystem which will provide 
%              momentum dumping capability for the spacecraft [kg] 
% F         = Thrust required to dump momentum 
% Hydrazine (monopropellant) is chosen as the fuel for this propulsion 
% system and a conservative specific impulse, Isp, is 200 seconds 
Isp = 200; 
% Here the earth's gravity constant is initialized (9.8 m/s^2) 
g = 9.8; 
% Here the impulse time, t, of the thruster firing is set 
% It is assumed that the thruster required for momentum dumping will fire 
% for 1 second 
t = 1; 
% Here the locations of the six required thrusters are initialized [x y z] 
% Each row is for a different thruster 
thruster = [0 size(2)/2 size(3)/2; 0 size(2)/2 -size(3)/2; size(1)/2 0 size(3)/2; -size(1)/2 ...
    0 size(3)/2; size(1)/2 -size(2)/2 0; -size(1)/2 -size(2)/2 0]; 
% Here the moment arms for the six thrusters from the CG are determined 
% For X-thrusters (spin about X-axis), moment arm is in Y-direction (cols 1,2) 
% For Y-thrusters (spin about Y-axis), moment arm is in Z-direction (cols 3,4) 
% For Z-thrusters (spin about Z-axis), moment arm is in X-direction (cols 5,6) 
moment_arms = [abs(CG(2) - thruster(1,2)) abs(CG(2)- thruster(2,2)) abs(CG(3) - ...
    thruster(3,3)) abs(CG(3) - thruster(4,3)) abs(CG(1) - thruster(5,1)) abs(CG(1) - thruster(6,1))]; 
% Here we will assume the worst-case distance from the thruster to the CG 
% (shortest) which will require the largest thrust to impart the required 
% torque on the S/C for momentum dumping 
worst_moment_arm = min(moment_arms); 
% Here the thrust required to dump the momentum is calculated (per pulse) 
F = H / (worst_moment_arm * t); 
% Here the required propellant mass for this propulsion system is estimated 
total_pulses = (lifetime * 365.25) / sat_rate; % total thruster pulses required over lifetime 
m_prop = (F * total_pulses * t)/(Isp * g);  % mass in kg 
% Here the total propulsion system mass is determined by assuming that 85% 
% of the propulsion system mass is propellant (SMAD, p. 660) - conservative 
p_mass = m_prop / 0.85;  % mass in kg