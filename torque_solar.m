% torque_solar.m 
function T_solar = torque_solar(A, CG, i, mat) 
% function to computer solar radiation pressure 
% INPUTS 
% A: vector describing size of object  
% CG: distance from center of solar pressure to center of mass (m) 
% i: angle of incidence of the Sun (radians) 
% mat: ID of material on outside of craft 
% OUTPUT 
% T_solar: solar radiation pressure, in N*m 
% some constants: 
% speed of light, m/s 
c = 3*10^8; 
% solar constant, W/m^2 
F_s = 1367;  
% get reflectance, q, from file based on material used 
tmp = readmatrix("C:\Users\sdbri\Downloads\Material_Properties_Table.csv");
q = tmp(mat,3); 

% find surface area of largest face of orbit 
A_s = A(1)*A(2); 
if(A(1)*A(3) > A_s) 
    A_s = A(1)*A(2); 
end 
if(A(2)*A(3) > A_s) 
    A_s = A(2)*A(3); 
end 
F = (F_s/c)*A_s*(1 + q)*cos(i); 
T_solar = F*(max(abs(CG))); 