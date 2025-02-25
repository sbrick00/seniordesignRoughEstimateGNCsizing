% oe2rv.m 
%  CREDIT: Christopher D. Hall  
%          http://www.aoe.vt.edu/~cdhall/ 
% 
%  oe2rv.m  Orbital Elements to r,v 
% 
%  [r,v] = oe2rv(oe,mu) 
%   
% oe = [a e i Om om nu] 
%   
% 
% r,v  expressed in  IJK  frame 
% a = semi-major axis 
% e = eccentricity 
% i = inclination 
% Om = argument of periapsis 
% om = right ascension of the ascending node (longitude of ascending node) 
% nu = true anomaly (at epoch). ***(location on orbit)*** 
function [ri,vi] = oe2rv(oe,mu) 
a=oe(1); e=oe(2); i=oe(3); Om=oe(4); om=oe(5); nu=oe(6); 
p = a*(1-e*e); 
r = p/(1+e*cos(nu)); 
rv = [r*cos(nu); r*sin(nu); 0];   
vv = sqrt(mu/p)*[-sin(nu); e+cos(nu); 0]; 
% 
% now rotate 
% 
cO = cos(Om);  sO = sin(Om); 
co = cos(om);  so = sin(om); 
ci = cos(i);   si = sin(i); 
R  = [cO*co-sO*so*ci  -cO*so-sO*co*ci  sO*si; 
  sO*co+cO*so*ci  -sO*so+cO*co*ci -cO*si; 
  so*si            co*si           ci]; 
ri = (R*rv)'; 
vi = (R*vv)';