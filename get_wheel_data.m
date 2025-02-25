% get_wheel_data.m 
function p = get_wheel_data 
wheel(1) = struct('name', 'Teldix RSI 01-5/15', 'ang_moment', 0.04, 'mass', 0.6); 
wheel(2) = struct('name', 'Teldix RSI 01-5/28', 'ang_moment', 0.12, 'mass', 0.7); 
wheel(3) = struct('name', 'LeoStar', 'ang_moment', 4.7, 'mass', 3.628); 
wheel(4) = struct('name', 'Dyncon MicroWheel 200', 'ang_moment', 0.18, 'mass', 0.93); 
wheel(5) = struct('name', 'Honeywell HR12', 'ang_moment', 50, 'mass', 9.5); 
wheel(6) = struct('name', 'Honeywell HR14', 'ang_moment', 75, 'mass', 10.6); 
wheel(7) = struct('name', 'Honeywell HR16', 'ang_moment', 100, 'mass', 12); 
wheel(8) = struct('name', 'Honeywell Miniature Reaction Wheel', 'ang_moment', 1.0, 'mass', 1.3); 
wheel(9) = struct('name', 'Honeywell HR0610', 'ang_moment', 12, 'mass', 5.0); 
wheel(10) = struct('name', 'Teldix DR23-0', 'ang_moment', 23, 'mass', 6.9); 
wheel(11) = struct('name', 'Teldix RDR68-6', 'ang_moment', 68, 'mass', 9.1); 
wheel(12) = struct('name', 'Teldix RSI 25-75/60', 'ang_moment', 25, 'mass', 6.3); 
wheel(13) = struct('name', 'Teldix RSI 68-75/60x', 'ang_moment', 68, 'mass', 8.5); 
wheel(14) = struct('name', 'Teldix RSI 4-75/60', 'ang_moment', 4, 'mass', 3.7); 
wheel(15) = struct('name', 'Teldix RSI 12-75/60x', 'ang_moment', 12, 'mass', 4.85); 
wheel(16) = struct('name', 'Teldix RSI 18-220/45', 'ang_moment', 18, 'mass', 6.3); 
wheel(17) = struct('name', 'Teldix RSI 30-280/30', 'ang_moment', 30, 'mass', 8.5); 
wheel(18) = struct('name', 'Teldix RSI 68-170/60', 'ang_moment', 68, 'mass', 8.9); 
wheel(19) = struct('name', 'Teldix RSI 02-25/30', 'ang_moment', 0.2, 'mass', 1.7); 
wheel(20) = struct('name', 'Teldix RSI 04-25/60', 'ang_moment', 0.4, 'mass', 1.7); 
wheel(21) = struct('name', 'Teldix RSI 1.6-25/60', 'ang_moment', 1.6, 'mass', 2.4); 
for(i=1:length(wheel)) 
    %plot(wheel(i).ang_moment, wheel(i).mass, 'r*'); 
    %hold on; 
    ang(i) = wheel(i).ang_moment; 
    mass(i) = wheel(i).mass; 
end 
[p,s] = polyfit(ang, mass, 4); 
%f = polyval(p, ang); 
%plot(ang, f, 'g*');