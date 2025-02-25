
%% Bill Nadir Disturbance Torques
%# Vehicle Properties #%
clc
close all
clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT INPUTS!!!!!!!!! - Ripped Straight from Bill Nadir. Use your
% Satellite values here. Check the PDF for info on material selection
% stuff. I did not use any of the reaction wheel or thruster stuff within
% this script tbh, but it's in there (Read the pdf lol).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

veh.dim = [1 1 1]; %[m], Length of edges on rectangular-prism shaped vehicle (length,
% width, depth)
veh.CG = [1 1 1]; %[m], center of gravity offset from geometric center
veh.mass = 718; %[kg] (minus ACS system)
veh.mat = 10; %Surface material code
veh.life = 5; %[yrs], Vehicle lifespan
%# Orbital Elements #%
OE.a = (500+6378)*1e3; %[m], semi-major axis
OE.e = 0.0; % eccentricity
OE.i = 97 *pi/180; %[rad], inclination
OE.Om = 0 *pi/180; %[rad], argument of periapsis (angle from ascending node to periapsis)
OE.om = 0; % [rad],longitude of the ascending node (angle between x and asc. node)
%# Planet Properties #%
planet.mu = 3.986e14; %[m^3/s^2], Earth gravity constant
planet.r_pol = 6357000; %[m], Polar radius
planet.r_equ = 6378000; %[m], Equitorial radius


%% Disturbance Torque Estimation Script - Ripped Straight From Bill Nadir
% I did not change anything about this. 

%# Calculate time for one orbit
t = 2*pi*sqrt(OE.a^3/planet.mu); %[sec]
%# Calculate Solar radiation torques
Ts = torque_solar(veh.dim, veh.CG, 0, veh.mat);
ang_step = pi/500; %[rad]
ang_range = 0:ang_step:2*pi;
ii = 1;
%## Calculate max disturbance torque around one complete orbit ##%
for ang = ang_range
    %# Calculate orbital postion and velocity
    [r,v] = oe2rv([OE.a OE.e OE.i OE.Om OE.om ang], planet.mu);
    v_mag = norm(v); %[m/s], Calculate velocity magnitude
    R(ii,:) = r; %[m]
    V(ii,:) = v; %[m/s]
    V_MAG(ii,:) = v_mag;
    E = acos((OE.e + cos(ang))/(1+OE.e*cos(ang))); %Eccentric anomoly
    time(ii,:) = sqrt(OE.a^3/planet.mu)*(E-OE.e*sin(E)); %time to 'ang'
    %# Post-process orbit elevation (latitude)
    aa = sqrt(r(1)^2 + r(2)^2);
    lat = atan2(r(3), aa); %[rad]
    LAT(ii,:) = lat;
    %# Post-process altitude
    r_planet = (planet.r_pol*planet.r_equ)/...
        sqrt(planet.r_pol^2*cos(lat)^2 + planet.r_equ^2*sin(lat)^2); %[m],

    alt = norm(r)-r_planet; %[m], Subtract planet radius from vehicle position vector
    % magnitude
    ALT(ii,:) = alt;
    %# Calculate Aerodynamic torques


    Ta = torque_aero(veh.dim, v_mag, alt, veh.CG);
    TA(ii,:) = Ta;
    %# Calculate Magnetic torques
    Tm = torque_magnetic(lat, norm(r), r_planet);
    TM(ii,:) = Tm;
    %# Calculate Gravity torques
    Tg = torque_gravity(norm(r), planet.mu, I);
    TG(ii,:) = Tg;
    TS(ii,:) = Ts;
    %# Sum all disturbance torques
    Ttot(ii,:) = Ts + Ta + Tm + Tg; %[Nm]
    ii=ii+1; %Increment counter
end
%## post-process time values ##%
max_t = ceil(length(time)/2);
for jj = [max_t+1: length(time)];
    time(jj) = 2*time(max_t)-time(jj);
end



%## Integrate max torques around orbit to find total ang mom ##% 
  ang_mom = trapz(time,Ttot); %[Nms], total angular momentum around one complete orbit 
  ang_mom_cyc = 0.8 * ang_mom; %[Nms], cyclical angular momentum per orbit 
  ang_mom_sec = ang_mom - ang_mom_cyc; %[Nms], Secular angular momentum per orbit 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % I DID NOT USE THIS SECTION! READ THE BILL NADIR PAPER FOR MORE INFO
  % I DID NOT USE THIS SECTION! READ THE BILL NADIR PAPER FOR MORE INFO
  % I DID NOT USE THIS SECTION! READ THE BILL NADIR PAPER FOR MORE INFO
  % I DID NOT USE THIS SECTION! READ THE BILL NADIR PAPER FOR MORE INFO


%## Size ACS actuators for cyclical momentum storage##% 
  wheel_data = get_wheel_data; %Loads Reaction wheel data (mass vs Nms)  
  w_mass = polyval(wheel_data, ang_mom_cyc); 
%## Size ACS thrusters for secular momentum dumping ##% 
   orb_sat = 1; %[orbits/saturation] 
   day_sat = orb_sat*t/86400; %[day/saturation] 
   [thrust, t_mass] = prop_system(veh.dim, veh.CG, veh.life, ang_mom_sec, day_sat); 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



   %Plotting 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % I USED THESE PLOTS IN MULTIPLE PRESENTATIONS. USE YOUR BEST JUDGEMENT 
   % I USED THESE PLOTS IN MULTIPLE PRESENTATIONS. USE YOUR BEST JUDGEMENT
   % I USED THESE PLOTS IN MULTIPLE PRESENTATIONS. USE YOUR BEST JUDGEMENT
   % I USED THESE PLOTS IN MULTIPLE PRESENTATIONS. USE YOUR BEST JUDGEMENT

%## Plot results ##% 
 % plot_planet_3D(R,planet); 
 figure(23); 
polarplot(ang_range', Ttot, 'k', LineWidth=1.15) 
 hold on 
 polarplot(ang_range', TG, LineWidth=1.15); 
polarplot(ang_range', TM,LineWidth=1.15); 
polarplot(ang_range', TA, LineWidth=1.1); 

polarplot(ang_range', TS, LineWidth=1.1); 
 hold off 
ggg = legend('Total', 'Gravity Gradient', 'Magnetic', 'Aero', 'Solar', 'FontSize', 12); 

% xlabel('Orbit Position From Ascending Node[deg]')

% Define custom labels while keeping orbital points labeled
theta_labels = string(0:30:360) + '°'; % Default labels
theta_labels([1, 7]) = {'Ascending Node', 'Descending Node'};
thetaticklabels(theta_labels);
pl = gca;
pl.ThetaAxis.FontSize = 12;
rticklabels({'0','1e-4 Nm', '2e-4 Nm', '3e-4 Nm'});
title('Extended Worst-Case Disturbance Torques', FontSize=14);
ax = gca; 
ax.RAxis.FontSize = 12;  % Change the radial label font size
ax.RAxis.FontWeight = "bold";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%% SIZING FOR SLEW REQUIREMENTS - Chat helped with this. 
% Separate from Bill Nadir script above
% I would upload this to chat and ask what it all means lol



% USER INPUTS
% J: 3x3 inertia matrix [kg m^2], assumed defined externally in your workspace
%    e.g., J = diag([Ixx, Iyy, Izz]) or a general 3x3. 
% Make sure it’s available before running this script.
if ~exist('J','var')
    error('Inertia matrix J is not defined in the workspace.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR SLEW REQUIREMENT INPUTS 
% YOUR SLEW REQUIREMENT INPUTS 
% YOUR SLEW REQUIREMENT INPUTS 
% YOUR SLEW REQUIREMENT INPUTS 
% YOUR SLEW REQUIREMENT INPUTS 
% YOUR SLEW REQUIREMENT INPUTS 

% Slew angle in degrees (example: 50 degrees)
slewAngle_deg = 15;  

% Slew time in seconds (example: 25 seconds)
slewTime      = 240;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% If you want to pick a specific principal axis of J (e.g., Z-axis), do:
%   axisOfInterest = 3;  % (1 for X, 2 for Y, 3 for Z)
% Or find the largest principal axis:
% [principalInertias, ~] = eig(J);   % For simple diagonal J, eigenvalues match the diag
% onal
% I_principal = max(diag(principalInertias));  % pick the largest principal moment
I_principal = max(Jbus, [], "all");  % pick the largest principal moment

% Or just pick the diagonal entry if J is diagonal and you know which axis you want:
% I_principal = J(3,3);  % for Z-axis, for instance.

% Convert slew angle from degrees to radians
slewAngle_rad = deg2rad(slewAngle_deg);

% Time split between accel & decel (bang-bang)
t_half = slewTime / 2;

% Angular acceleration for a symmetric trapezoid:
%   total angle = alpha * t_half^2 (accel) + alpha * t_half^2 (decel) = alpha * t_half^2 * 2
%   => slewAngle_rad = alpha * t_half^2
alpha = slewAngle_rad / (t_half^2);   % [rad/s^2]

% Peak angular velocity at midpoint
omega_max = alpha * t_half;           % [rad/s]

% Required torque about that axis
tau_required = I_principal * alpha;   % [N·m]

% Peak spacecraft momentum about that axis
H_spacecraft = I_principal * omega_max; % [N·m·s]

% 3) REPORT THE MINIMUM NEEDED SPECS

% For a Reaction Wheel, you need:
% - Enough torque to accelerate the spacecraft at 'tau_required' (though many RW specs
%   center on momentum rather than continuous torque).
% - Enough momentum capacity to store 'H_spacecraft' during the slew.

min_RW_momentum_needed = H_spacecraft;    % [N·m·s]
min_RW_torque_needed   = tau_required;    % [N·m] (if you want a direct torque rating)

% For a CMG (single-gimbal or array), you need:
% - Max torque >= tau_required
% - Rotor momentum capacity >= H_spacecraft

min_CMG_momentum_needed = H_spacecraft;   % [N·m·s]
min_CMG_torque_needed   = tau_required;   % [N·m]

% 4) DISPLAY RESULTS

fprintf('=== Maneuver Requirements ===\n');
fprintf('Slew Angle       : %.2f deg (%.4f rad)\n', slewAngle_deg, slewAngle_rad);
fprintf('Slew Time        : %.2f s\n', slewTime);
fprintf('Principal Inertia: %.3f kg·m^2\n', I_principal);
fprintf('Angular Accel    : %.4e rad/s^2\n', alpha);
fprintf('Peak Angular Rate: %.4e rad/s (%.4f deg/s)\n', omega_max, rad2deg(omega_max));
fprintf('Required Torque  : %.4e N·m\n', tau_required);
fprintf('Peak Momentum    : %.4e N·m·s\n\n', H_spacecraft);

fprintf('=== Minimum Specifications Needed ===\n');
fprintf('Reaction Wheel:\n');
fprintf('  - Momentum Capacity >= %.3f N·m·s\n', min_RW_momentum_needed);
fprintf('  - Torque Rating    >= %.3f N·m (optional check)\n\n', min_RW_torque_needed);

fprintf('CMG:\n');
fprintf('  - Momentum Capacity >= %.3f N·m·s\n', min_CMG_momentum_needed);
fprintf('  - Max Torque        >= %.3f N·m\n', min_CMG_torque_needed);

fprintf('\nNote: Add margins for disturbances, uncertainties, and design safety.\n');
