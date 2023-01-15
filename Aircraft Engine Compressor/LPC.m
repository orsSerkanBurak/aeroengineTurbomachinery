clc;clear;close all;

%% Aerodynamic Definition Calculations of LPC

D=0.5; % Diffusion factor
M1 = 0.4; % Fan input velocity in Mach # 
gamma = 1.4; % Ratio of specific heats
sigma = 1; % Solidity of fan
ec = 0.9039; % Polytropic efficiency of fan
alpha1 = -89.9:0.1:89.9; % blade angle
Gamma = ((2*sigma)+(sind(alpha1)))./(cosd(alpha1)); 
alpha2 = acosd(((2*sigma.*Gamma*(1-D))+(sqrt((Gamma.^2)+1-((4*sigma^2)*(1-D)^2))))./((Gamma.^2)+1));
deltaAlpha = alpha2 - alpha1; % Turning angle
wr_V1=(cosd(alpha1)).*((tand(alpha1))+((tand(alpha2)))); % Rotational speed-inlet velocity ratio
M1R_M1=(cosd(alpha1))./(cosd(alpha2)); % relative velocity - inlet velocity ratio
tauS = ((((gamma-1)*M1^2)/(1+((gamma-1)*((M1^2)/2))))*(((cosd(alpha1)).^2./(cosd(alpha2)).^2)-1))+1; % Stage temperature ratio
piS = tauS.^((gamma*ec)/(gamma-1)); % Stage pressure ratio
deltaTt = (((cos(alpha1)).^2)./((cos(alpha2)).^2))-1; % Stage total temperature increase
Tt1 = deltaTt./(tauS-1);
M1R = M1R_M1.*M1; % Inlet relative Mach number

%blade angles
dc = (alpha1-alpha2)./((4*sqrt(sigma))-1);
gamma1 = alpha1;
gamma2 = alpha2 - dc;

%plot
figure(1)
yyaxis left
p1 = plot(alpha1,wr_V1,'Color',	[0, 0.75, 0.75],'LineStyle','-','LineWidth',1.2);
hold on;
p2 = plot(alpha1,M1R_M1,'Color',[0.75, 0.75, 0],'LineStyle','-','LineWidth',1.2);
p3 = plot(alpha1,piS,'Color',[0.4940, 0.1840, 0.5560],'LineStyle','-','LineWidth',1.2);
xlabel('\alpha_{1} (deg)');
ylabel('\omegar / V_{1}      \pi_{s}      M_{1R} / M_{1}');
axis([-90 90 -1 4]);
yyaxis right 
p4 = plot(alpha1,alpha2,'Color',[0.6350, 0.0780, 0.1840],'LineStyle','-','LineWidth',1.2);
p5 = plot(alpha1,deltaAlpha,'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-','LineWidth',1.2);
ylabel('\alpha_{2} (deg)               \Delta\alpha (deg)');
grid on;
title('Repeating Compressor Stage(D=0.5,M_{1}=0.4,\gamma =1.4,\sigma=1,e_{c} = 0.9039)');
legend('\omegar / V_{1}','M_{1R} / M_{1}','\pi_{s}','\alpha_{2}','\Delta\alpha','Location','eastoutside');
axis([-90 90 0 90]);
hold off;
