clc;clear;close all;
%% HPT Aerodynamic Calculations

alpha2 = -80:0.1:80; % Turning angle
M3R = 0.9; % Relative Mach number
M2 = 1.1; % Stage exit Mach number
gamma_t = 1.3; % ratio of specific heats
omega = 0.3449; % dimensionless turbine rotor speed

% Velocity Ratios
V2overV = sqrt(((gamma_t-1)*M2.^2)/(1+((gamma_t-1)*M2.^2)/2));
uOverV = V2overV*cosd(alpha2);
v2OverV = V2overV*sind(alpha2);
v2ROverV = v2OverV - omega;

Beta2 = atand(v2ROverV./uOverV);
rotorRelative = 1 + ((omega^2)*(0.5-(v2OverV./omega)));
Beta3 = atand((rotorRelative./uOverV.^2)*(((gamma_t-1)*M3R^2)/(1+((gamma_t-1)*M3R^2)/2)));
Beta = Beta2 + Beta3;
v3OverV = (uOverV.*tand(Beta3))-omega;
alpha3 = atand(v3OverV./uOverV);
T3OverTt1 = (rotorRelative)/(1+(((gamma_t-1)*M3R^2)/2)); % Temperature ratio
tauts = T3OverTt1 + ((uOverV.^2).*((1+(tand(alpha3)).^2)/2)); % Stage total temperature ratio
solidity = 2*((cosd(Beta3)).^2).*((tand(Beta2))+((tand(Beta3)))); % Rotor solidity

% Degree of Reaction
Tt1 = 3162.29;
gcOverCpt = 7378;
V = sqrt(gcOverCpt*Tt1);
V3R = 1004.8;
v3r = V3R*sin(Beta3);
v2r = (v2OverV-omega)*V;
Rt = ((v3r.^2)-(v2r.^2))./(2*(1-tauts)*V^2); % Rotor degree of reaction


figure(1)
plot(alpha2,tauts,'LineWidth',1.2);
xlabel('\alpha_{2} (deg)'); ylabel('\tau_{ts}');
grid on;
title('Stage total temperature ratio at \Omega = 0.3449, M_{2} = 1.1, M_{3R} = 0.9 and \gamma_{t} = 1.3');

figure(2)
plot(alpha2,Beta,'LineWidth',1.2);
xlabel('\alpha_{2} (deg)'); ylabel('\beta_{2}+\beta_{3} (deg)');
grid on;
title('Rotor flow turning angle at \Omega = 0.3449, M_{2} = 1.1, M_{3R} = 0.9 and \gamma_{t} = 1.3');

figure(3)
plot(alpha2,alpha3,'LineWidth',1.2);
xlabel('\alpha_{2} (deg)'); ylabel('\alpha_{3} (deg)');
grid on;
title('Stage exit flow angle at \Omega = 0.3449, M_{2} = 1.1, M_{3R} = 0.9 and \gamma_{t} = 1.3');

figure(4)
plot(alpha2,solidity,'LineWidth',1.2);
xlabel('\alpha_{2} (deg)'); ylabel('\sigma_{xr}');
grid on;
title('Rotor solidity based on Z=1 at \Omega = 0.3449, M_{2} = 1.1, M_{3R} = 0.9 and \gamma_{t} = 1.3');

figure(5)
plot(alpha2,Rt,'LineWidth',1.2);
xlabel('\alpha_{2} (deg)'); ylabel('R_{t}');
grid on;
title('Rotor degree of reaction at \Omega = 0.3449, M_{2} = 1.1, M_{3R} = 0.9 and \gamma_{t} = 1.3');
axis([-80 80 -1 1.5]);