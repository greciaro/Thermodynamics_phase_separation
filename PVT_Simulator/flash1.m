function [xi,yi,Ki,VL,VV,l,v]=FLASH_A_Astro(P,T,zi,Tci,Pci,wi) 
% INPUT: Total composition.
% OUTPUT: Liquid and Vapor composition and equilibrium constants.

% clear all

% 1.- Calculating initial values and program parameters.
R = 8.31438;%Pa-m^3/mol-K
% P = 2000;%psia
% P = P*6892.8571;%Pa
% T = 162;%F
% T = (T-32)/1.8+273.15;%K
% Tci = [-116.63 305.65 652.1 87.9];%F
% Tci = ((Tci-32)./1.8)+273.15;%K
% Pci = [667.8 550.7 305.7 1071];%psia
% Pci = Pci.*6892.8571;%Pa
% wi = [.0104 .201 .49 .225];%Ascentric factor.
%--------------------------------------------------------------
% T = 162;%F
% T = (T-32)/1.8+273.15;%K
% P = 4.24*10^6;
% Tci = zeros(1,12);
% Pci = zeros(1,12);
% wi = zeros(1,12);
% zi = zeros(1,12);
% 
% Mole_fraction = [1.150 0.208 34.735 11.072 9.169 5.508 3.405 4.490 4.560 3.720 2.540 2.265];%percentage
% zi = (Mole_fraction./sum(Mole_fraction));%fraction
% 
% % 1.- N2
% Tci(1) = 126.2; 
% Pci(1) = 33.9*10^5; 
% wi(1) = 0.039; 
% % 2.- CO2
% Tci(2) = 304.1;
% Pci(2) = 73.8*10^5;
% wi(2) = 0.239;
% % 3.- C1
% Tci(3) = 190.4;
% Pci(3) = 46*10^5;
% wi(3) = 0.011;
% % 4.- C2
% Tci(4) = 305.4;
% Pci(4) = 48.8*10^5;
% wi(4) = 0.099;
% % 5.- C3
% Tci(5) = 369.8;
% Pci(5) = 42.5*10^5;
% wi(5) = 0.153;
% % 6.- C4
% Tci(6) = 425.2;
% Pci(6) = 38*10^5;
% wi(6) = 0.199;%nC4
% % 7.- C5
% Tci(7) = 469.7;
% Pci(7) = 33.7*10^5;
% wi(7) = 0.251;%nC5
% % 8.- C6
% Tci(8) = 507.5;
% Pci(8) = 30.1*10^5;
% wi(8) = 0.299; 
% % 9.- C7
% Tci(9) = 540.3;
% Pci(9) = 27.4*10^5;
% wi(9) = 0.349; 
% % 10.- C8
% Tci(10) = 568.8;
% Pci(10) = 24.9*10^5;
% wi(10) = 0.398;
% % 11.- C9
% Tci(11) = 594.6;
% Pci(11) = 22.9*10^5;
% wi(11) = 0.445;
% % 12.- C10
% Tci(12) = 617.7;
% Pci(12) = 21.2*10^5;
% wi(12) = 0.489;

%-----------------------------------------------------------------

Tr = T./Tci;%Reduced temperatue
E = ones(1,length(Tci))*1e-6;%Tolerance to reach Ki for equilibrium.
dKi = ones(1,length(Tci))*1000;%intial value to compare tolerance.
absdKi = abs(dKi);%intial value to compare tolerance.
Ki = (Pci./P).*exp(5.37.*(1+wi).*(1-(Tci./T)));%initial value of Ki with Wilson equation.
alfa = 0.5;%proportion to mix liquid and vapor and get zi.
%%
