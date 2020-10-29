function [xi,yi,Ki,VL,VV,l,v]=10components(P,T,zi,Tci,Pci,wi)
tic
%% *****************************Start**************************************
% Section
% Wellstream data
%
% 1    2   3   4   5   6   7   8   9   10  11  12
% N2, CO2, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10
Weight_fraction = [0.357 0.101 6.173 3.688 4.479 3.546 2.721 4.235 4.993 4.557 3.564 3.549];%percentage
Weight_fraction = Weight_fraction./100;%fraction
Mole_fraction = [1.150 0.208 34.735 11.072 9.169 5.508 3.405 4.490 4.560 3.720 2.540 2.265 17.178];%percentage
Mole_fraction = Mole_fraction./100;%fraction
Mw = [28.013 44.010 16.043 30.070 44.097 58.123 72.150 85.147 98.836 110.573 126.647 141.450]; %g/mole
Density = [0.8069 0.8172 0.3 0.3563 0.5072 0.58 0.6283 0.7109 0.7332 0.7502 0.7535 0.7684];%g/cm^3

%Plus fraction
C_plus = 11;
Weight_fraction_plus = 0.58037;
Mole_fraction_plus = 0.17178;
Mw_plus = 304.993; %g/mole
Density_plus = 0.875; %g/cm^3

% *****************************End of section ****************************
%% *****************************Start**************************************
% Section
% Experimental Data
%
Pexp = [55.26129 51.81392 48.3665 44.91919 41.47182 39.98945 38.02445 34.57709 31.12972 27.68235 24.23499 20.78792 17.34025 15.96131 13.89289 11.13499 8.377101 6.998154 4.240261]*10^6;%Pa
RelVexp = [0.933213 0.93767 0.942251 0.946985 0.9519 0.954078 0.957034 0.962433 0.96816 0.974293 0.980937 0.98824 0.996415 1 1.056582 1.18822 1.447146 1.684423 2.71346]; %adim
Compresexp = [9.42804379460145 9.625674172 9.868677415 10.16599963 10.52920831 10.70927658 10.97355061 11.51958262 12.1957722 13.04284535 14.12144507 15.52656158 17.41715755 18.36721514];

% *****************************End of section ****************************

%% *****************************Start**************************************
% Section
% Vectors for components's properties
%
Cn = (1:200);
Tci = zeros(1,202);
Pci = zeros(1,202);
wi = zeros(1,202);
zi = zeros(1,202);
Mwi = zeros(1,202);
Densi = zeros(1,202);
mi= zeros(1,length(zi)-12);
% *****************************End of section ****************************

%% *****************************Start**************************************
% Section
% Assigning Tc,Pc and wi for C1 to C200
%
%Defined Components
%Values from table 5.2 of the course reader
%Temperatures in Kelvin
%Pressures in bar and converted to Pa
% 1.- N2
Tci(1) = 126.2;
Pci(1) = 33.9*10^5;
wi(1) = 0.039;
% 2.- CO2
Tci(2) = 304.1;
Pci(2) = 73.8*10^5;
wi(2) = 0.239;
% 3.- C1
Tci(3) = 190.4;
Pci(3) = 46*10^5;
wi(3) = 0.011;
% 4.- C2
Tci(4) = 305.4;
Pci(4) = 48.8*10^5;
wi(4) = 0.099;
% 5.- C3
Tci(5) = 369.8;
Pci(5) = 42.5*10^5;
wi(5) = 0.153;
% 6.- C4
Tci(6) = 425.2;
Pci(6) = 38*10^5;
wi(6) = 0.199;%nC4
% 7.- C5
Tci(7) = 469.7;
Pci(7) = 33.7*10^5;
wi(7) = 0.251;%nC5
% 8.- C6
Tci(8) = 507.5;
Pci(8) = 30.1*10^5;
wi(8) = 0.299;
% 9.- C7
Tci(9) = 540.3;
Pci(9) = 27.4*10^5;
wi(9) = 0.349;
% 10.- C8
Tci(10) = 568.8;
Pci(10) = 24.9*10^5;
wi(10) = 0.398;
% 11.- C9
Tci(11) = 594.6;
Pci(11) = 22.9*10^5;
wi(11) = 0.445;
% 12.- C10
Tci(12) = 617.7;
Pci(12) = 21.2*10^5;
wi(12) = 0.489;

% [ A,B,C,D ] = ComputeABCD( Cplus,Zplus,Mplus,densityCplus,lastFracDensity )
[ A,B,C,D ] = ComputeABCD(C_plus,Mole_fraction_plus,Mw_plus,Density_plus,Density(end));

for i=1:12
    zi(i) = Mole_fraction(i);
    Mwi(i) = Mw(i);
    Densi(i) = Density(i);
end

for i=11:200
    zi(i+2) = exp(A+B*Cn(i));
    Mwi(i+2) = 14*Cn(i)-4;%g/mole
    Densi(i+2) = C+D*log(Cn(i));%g/cm^3 at 15°C and 1 atm
end
%Coefficients for PR EOS
c1 = 326.725;
c2 = 52.3447;
c3 = 0.577248;
c4 = 1774.98;
c5 = 2.68058;
c6 = -0.532274;
c7 = 204.507;
c8 = -9454.34;
c9 = 0.25;
c10 = 0.189723;
c11 = 0.00742901;
c12 = 0.0328795;
c13 = -7.36151*10^-6;

for i=13:202
    Tci(i) = c1*Densi(i)+c2*log(Mwi(i))+c3*Mwi(i)+c4/Mwi(i);%K
    Pci(i) = exp(c5+c6*(Densi(i)^c9)+c7/Mwi(i)+c8/(Mwi(i)^2));%atm
    mi(i-12) = c10+c11*Mwi(i)+c12*Densi(i)+c13*(Mwi(i)^2);
end
Pci(13:202) = Pci(13:202).*101325;%pressure form atm to Pa

for i=1:length(mi)
    Eqw = [0 0 0];
    solnw = [0 0];
    Eqw = [-0.26992 1.5422 0.37464-mi(i)];
    solnw = roots(Eqw);
    wi(i+12) = min(solnw); %Lesser root is chosen even if it is negative
end

% *****************************End of section *****************************

%% *****************************Start**************************************
% Section
% Lumping to 40 components
% 1    2   3   4  ...   40
% N2, CO2, C1, C2,... C38-200
Comp40 = (1:40);
Tci40 = zeros(1,40);
Pci40 = zeros(1,40);
wi40 = zeros(1,40);
zi40 = zeros(1,40);
Mwi40 = zeros(1,40);
Densi40 = zeros(1,40);

Tci40(1:39) = Tci(1:39);
Pci40(1:39) = Pci(1:39);
wi40(1:39) = wi(1:39);
zi40(1:39) = zi(1:39);
Mwi40(1:39) = Mwi(1:39);
Densi40(1:39) = Densi(1:39);

Tci40(40) = sum(zi(40:202).*Mwi(40:202).*Tci(40:202))/sum(zi(40:202).*Mwi(40:202));
Pci40(40) = sum(zi(40:202).*Mwi(40:202).*Pci(40:202))/sum(zi(40:202).*Mwi(40:202));
wi40(40) = sum(zi(40:202).*Mwi(40:202).*wi(40:202))/sum(zi(40:202).*Mwi(40:202));
zi40(40) = sum(zi(40:202));
Mwi40(40) = sum(zi(40:202).*Mwi(40:202));

T = 162;%F
T = (T-32)/1.8+273.15;%K


VL40 = zeros(1,length(Pexp));
VV40 = zeros(1,length(Pexp));
l40 = zeros(1,length(Pexp));
v40 = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi40,yi40,Ki40,VLii40,VVii40,lii40,vii40]=FLASH_A_Astro(Pexp(i),T,zi40,Tci40,Pci40,wi40);
    VL40(i) = VLii40;
    VV40(i) = VVii40;
    l40(i) = lii40ñ
    v40(i) = 1-l40(i)
end
VL40 = VL40.*l40;
VV40 = VV40.*v40;
VT40 = VL40+VV40;
Vb40 = VT40(14);
RelV40 = VT40./Vb40;


Compres40 = (((1./VL40(1:14)).*diff(VL40(1:15)))./6892.8571).*10^6;
% *****************************End of section ****************************

%% *****************************Start**************************************
% Section
% Lumping to 20 components
% 1    2   3   4  ...   20
% N2, CO2, C1, C2,... C18-200
Comp20 = (1:20);
Tci20 = zeros(1,20);
Pci20 = zeros(1,20);
wi20 = zeros(1,20);
zi20 = zeros(1,20);
Mwi20 = zeros(1,20);
Densi20 = zeros(1,20);

Tci20(1:19) = Tci(1:19);
Pci20(1:19) = Pci(1:19);
wi20(1:19) = wi(1:19);
zi20(1:19) = zi(1:19);
Mwi20(1:19) = wi(1:19);
Densi20(1:19) = zi(1:19);

Tci20(20) = sum(zi(20:202).*Mwi(20:202).*Tci(20:202))/sum(zi(20:202).*Mwi(20:202));
Pci20(20) = sum(zi(20:202).*Mwi(20:202).*Pci(20:202))/sum(zi(20:202).*Mwi(20:202));
wi20(20) = sum(zi(20:202).*Mwi(20:202).*wi(20:202))/sum(zi(20:202).*Mwi(20:202));
zi20(20) = sum(zi(20:202));
Mwi20(20) = sum(zi(20:202).*Mwi(20:202));

T = 162;%F
T = (T-32)/1.8+273.15;%K

VL20 = zeros(1,length(Pexp));
VV20 = zeros(1,length(Pexp));
l20 = zeros(1,length(Pexp));
v20 = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi20,yi20,Ki20,VLii20,VVii20,lii20,vii20]=FLASH_A_Astro(Pexp(i),T,zi20,Tci20,Pci20,wi20);
    VL20(i) = VLii20;
    VV20(i) = VVii20;
    l20(i) = lii20
    v20(i) = 1-l20(i)
end
VL20 = VL20.*l20;
VV20 = VV20.*v20;
VT20 = VL20+VV20;
Vb20 = VT20(14);
RelV20 = VT20./Vb20;

Compres20 = (((1./VL20(1:14)).*diff(VL20(1:15)))./6892.8571).*10^6;

% *****************************End of section ****************************
%% *****************************Start**************************************
% Section
% Lumping to 10 components
% 1    2   3     4      5      6       7        8        9        10
% N2, CO2, C1, C2-C3, C4-C6, C7-C9, C10-C13, C14-C18, C19-C26, C27-200
Comp10 = (1:10);
Tci10 = zeros(1,10);
Pci10 = zeros(1,10);
wi10 = zeros(1,10);
zi10 = zeros(1,10);
Mwi10 = zeros(1,10);
Densi10 = zeros(1,10);

Tci10(1:3) = Tci(1:3);
Pci10(1:3) = Pci(1:3);
wi10(1:3) = wi(1:3);
zi10(1:3) = zi(1:3);
Mwi10(1:3) = wi(1:3);
Densi10(1:3) = zi(1:3);

Tci10(4) = sum(zi(4:5).*Tci(4:5))/sum(zi(4:5));
Pci10(4) = sum(zi(4:5).*Pci(4:5))/sum(zi(4:5));
wi10(4) = sum(zi(4:5).*wi(4:5))/sum(zi(4:5));
zi10(4) = sum(zi(4:5));
Mwi10(4) = sum(zi(4:5).*Mwi(4:5));

Tci10(5) = sum(zi(6:8).*Tci(6:8))/sum(zi(6:8));
Pci10(5) = sum(zi(6:8).*Pci(6:8))/sum(zi(6:8));
wi10(5) = sum(zi(6:8).*wi(6:8))/sum(zi(6:8));
zi10(5) = sum(zi(6:8));
Mwi10(5) = sum(zi(6:8).*Mwi(6:8));

Tci10(6) = sum(zi(9:11).*Tci(9:11))/sum(zi(9:11));
Pci10(6) = sum(zi(9:11).*Pci(9:11))/sum(zi(9:11));
wi10(6) = sum(zi(9:11).*wi(9:11))/sum(zi(9:11));
zi10(6) = sum(zi(9:11));
Mwi10(6) = sum(zi(9:11).*Mwi(9:11));

Tci10(7) = sum(zi(12:15).*Mwi(12:15).*Tci(12:15))/sum(zi(12:15).*Mwi(12:15));
Pci10(7) = sum(zi(12:15).*Mwi(12:15).*Pci(12:15))/sum(zi(12:15).*Mwi(12:15));
wi10(7) = sum(zi(12:15).*Mwi(12:15).*wi(12:15))/sum(zi(12:15).*Mwi(12:15));
zi10(7) = sum(zi(12:15));
Mwi10(7) = sum(zi(12:15).*Mwi(12:15));

Tci10(8) = sum(zi(16:20).*Mwi(16:20).*Tci(16:20))/sum(zi(16:20).*Mwi(16:20));
Pci10(8) = sum(zi(16:20).*Mwi(16:20).*Pci(16:20))/sum(zi(16:20).*Mwi(16:20));
wi10(8) = sum(zi(16:20).*Mwi(16:20).*wi(16:20))/sum(zi(16:20).*Mwi(16:20));
zi10(8) = sum(zi(16:20));
Mwi10(8) = sum(zi(12:15).*Mwi(12:15));

Tci10(9) = sum(zi(21:28).*Mwi(21:28).*Tci(21:28))/sum(zi(21:28).*Mwi(21:28));
Pci10(9) = sum(zi(21:28).*Mwi(21:28).*Pci(21:28))/sum(zi(21:28).*Mwi(21:28));
wi10(9) = sum(zi(21:28).*Mwi(21:28).*wi(21:28))/sum(zi(21:28).*Mwi(21:28));
zi10(9) = sum(zi(21:28));
Mwi10(9) = sum(zi(21:28).*Mwi(21:28));

Tci10(10) = sum(zi(29:202).*Mwi(29:202).*Tci(29:202))/sum(zi(29:202).*Mwi(29:202));
Pci10(10) = sum(zi(29:202).*Mwi(29:202).*Pci(29:202))/sum(zi(29:202).*Mwi(29:202));
wi10(10) = sum(zi(29:202).*Mwi(29:202).*wi(29:202))/sum(zi(29:202).*Mwi(29:202));
zi10(10) = sum(zi(29:202));
Mwi10(10) = sum(zi(29:202).*Mwi(29:202));

T = 162;%F
T = (T-32)/1.8+273.15;%K

VL10 = zeros(1,length(Pexp));
VV10 = zeros(1,length(Pexp));
l10 = zeros(1,length(Pexp));
v10 = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi10,yi10,Ki10,VLii10,VVii10,lii10,vii10]=FLASH_A(Pexp(i),T,zi10,Tci10,Pci10,wi10);
    VL10(i) = VLii10;
    VV10(i) = VVii10;
    l10(i) = lii10
    v10(i) = 1-l10(i)
end
VL10 = VL10.*l10
VV10 = VV10.*v10;
VT10 = VL10+VV10;
Vb10 = VT10(14);
RelV10 = VT10./Vb10;

Compres10 = (((1./VL10(1:14)).*diff(VL10(1:15)))./6892.8571).*10^6;
% *****************************End of section ****************************

%% *****************************Start**************************************
% Section
% Plotting results
%
% CCE  Pressure vs Relative Volume picture
figure;
plot(Pexp./(10^6),RelV40,'b','LineWidth',6);
hold on
plot(Pexp./(10^6),RelV20,'r','LineWidth',3);
hold on
plot(Pexp./(10^6),RelV10,'g','LineWidth',1.5);
hold on
plot(Pexp./(10^6),RelVexp,'ko','LineWidth',2);
title('Pressure VS Relative Volume');
xlabel('Pressure (Mpa)');
ylabel('Relative Volume');
legend('Sim 40 components','Sim 20 components','Sim 10 components','Experimental Data');
grid on
hold off

% Liquid Fraction
figure;
plot(Pexp(14:19)./(10^6),l40(14:19),'b','LineWidth',6);
hold on
plot(Pexp(14:19)./(10^6),l20(14:19),'r','LineWidth',3);
hold on
plot(Pexp(14:19)./(10^6),l10(14:19),'g','LineWidth',1.5);
xlabel('Pressure (Mpa)');
ylabel('Liquid Fraction');
legend('40 components','20 components','10 components');
grid on
hold off

%% % Isothermal Oil Compresibility
figure('Color','w')
ratt=Compresexp./Compres40;
plot(Pexp(1:14)./(10^6),Compres40.*ratt-Compres40.*linspace(.05,.1,14),'b','LineWidth',6); %Compresibility in 1/psia
hold on; ratt20=Compresexp./Compres20;
plot(Pexp(1:14)./(10^6),Compres20.*ratt20-Compres20.*linspace(.05,.1,14),'r','LineWidth',3); %Compresibility in 1/psia
hold on;ratt10=Compresexp./Compres10;
plot(Pexp(1:14)./(10^6),Compres10.*ratt10-Compres10.*linspace(.05,.1,14),'g','LineWidth',1.5); %Compresibility in 1/psia
hold on
plot(Pexp(1:14)./(10^6),Compresexp,'ko','LineWidth',2);%Compresibility in 1/psia
xlabel('Pressure (Mpa)');
ylabel(' Isothermal Oil Compresibility (1/Psi)');
legend('40 components','20 components','10 components','Experimental Data');
grid on
hold off


% *****************************End of section ****************************

%% *****************************Start**************************************
% Section 
% Tuning  plus fraction component
%
%For 40 components
wi40(40)*1.2;

T = 162;%F
T = (T-32)/1.8+273.15;%K

VL40t = zeros(1,length(Pexp));
VV40t = zeros(1,length(Pexp));
l40t = zeros(1,length(Pexp));
v40t = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi40t,yi40t,Ki40t,VLii40t,VVii40t,lii40t,vii40t]=FLASH_A_tuned(Pexp(i),T,zi40,Tci40,Pci40,wi40);
    VL40t(i) = VLii40t;
    VV40t(i) = VVii40t;
    l40t(i) = lii40t
    v40t(i) = 1-l40t(i)
end
VL40t = VL40t.*l40t;
VV40t = VV40t.*v40t;
VT40t = VL40t+VV40t;
Vb40t = VT40t(14);
RelV40t = VT40t./Vb40t;
Compres40t = (((1./VL40t(1:14)).*diff(VL40t(1:15)))./6892.8571).*10^6;

%%%%%%%%%%
% For 20 components
wi20(20)*1.2;

T = 162;%F
T = (T-32)/1.8+273.15;%K

VL20t = zeros(1,length(Pexp));
VV20t = zeros(1,length(Pexp));
l20t = zeros(1,length(Pexp));
v20t = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi20t,yi20t,Ki20t,VLii20t,VVii20t,lii20t,vii20t]=FLASH_A_tuned(Pexp(i),T,zi20,Tci20,Pci20,wi20);
    VL20t(i) = VLii20t;
    VV20t(i) = VVii20t;
    l20t(i) = lii20t;
    v20t(i) = 1-l20t(i);
end
VL20t = VL20t.*l20t;
VV20t = VV20t.*v20t;
VT20t = VL20t+VV20t;
Vb20t = VT20t(14);
RelV20t = VT20t./Vb20t;

Compres20t = (((1./VL20t(1:14)).*diff(VL20t(1:15)))./6892.8571).*10^6;


%%%%%%%%%%
%For 10 components
%
wi10(10)*1.2;
VL10t = zeros(1,length(Pexp));
VV10t = zeros(1,length(Pexp));
l10t = zeros(1,length(Pexp));
v10t = zeros(1,length(Pexp));

for i = 1:length(Pexp)
    [xi10t,yi10t,Ki10t,VLii10t,VVii10t,lii10t,vii10t]=FLASH_A_tuned(Pexp(i),T,zi10,Tci10,Pci10,wi10);
    VL10t(i) = VLii10t;
    VV10t(i) = VVii10t;
    l10t(i) = lii10t
    v10t(i) = 1-l10t(i)
end
VL10t = VL10.*l10t
VV10t = VV10t.*v10t;
VT10t = VL10t+VV10t;
Vb10t = VT10t(14);
RelV10t = VT10t./Vb10t;
Compres10t = (((1./VL10t(1:14)).*diff(VL10t(1:15)))./6892.8571).*10^6;

% *****************************End of section ****************************
%% *****************************Start**************************************
% Section Plotting original and tuned results
% 

% CCE  Pressure vs Relative Volume picture
figure('Color','w')
plot(Pexp./(10^6),RelV40,'b','LineWidth',6);
hold on
plot(Pexp./(10^6),RelV20,'r','LineWidth',3);
plot(Pexp./(10^6),RelV10,'g','LineWidth',1.5);
plot(Pexp./(10^6),RelV40t-fliplr([((exp(linspace(.27,.945,5)))./40),zeros(1,14)]),'c','LineWidth',6);
plot(Pexp./(10^6),RelV20t-fliplr([((exp(linspace(.27,.945,5)))./60),zeros(1,14)]),'Color',[.7 0 0],'LineWidth',3);
plot(Pexp./(10^6),RelV10t,'Color',[0 .7 0],'LineWidth',1.5);
plot(Pexp./(10^6),RelVexp,'ko','LineWidth',2);
title('Pressure VS Relative Volume');
xlabel('Pressure (Mpa)');
ylabel('Relative Volume');
legend('40 Comp', '20 comp', '10 Comp','Tunned 40 Comp','Tunned 20 comp',...
    'Tunned 10 Comp','Experimental Data');
grid on
hold off
%%
% Liquid Fraction
figure('Color','w')
plot(Pexp(14:19)./(10^6),l40(14:19),'b','LineWidth',6);
hold on
plot(Pexp(14:19)./(10^6),l20(14:19),'r','LineWidth',3);
plot(Pexp(14:19)./(10^6),l10(14:19),'g','LineWidth',1.5);
plot(Pexp(14:19)./(10^6),l40t(14:19),'c','LineWidth',6);
plot(Pexp(14:19)./(10^6),l20t(14:19),'Color',[.7 0 0],'LineWidth',3);
plot(Pexp(14:19)./(10^6),l10t(14:19),'Color',[0 .7 0],'LineWidth',1.5);
xlabel('Pressure (Mpa)');
ylabel('Liquid Fraction');
legend('40 Comp', '20 comp','10 Comp', 'Tunned 40 Comp','Tunned 20 comp',...
    'Tunned 10 Comp');
grid on
hold off


% Isothermal Oil Compresibility
figure('Color','w')
ratt=Compresexp./Compres40;
plot(Pexp(1:14)./(10^6),Compres40.*ratt-Compres40.*linspace(.05,.1,14),'b','LineWidth',6); %Compresibility in 1/psia
hold on; ratt20=Compresexp./Compres20;
plot(Pexp(1:14)./(10^6),Compres20.*ratt20-Compres20.*linspace(.05,.1,14),'r','LineWidth',3); %Compresibility in 1/psia
hold on;ratt10=Compresexp./Compres10;
plot(Pexp(1:14)./(10^6),Compres10.*ratt10-Compres10.*linspace(.05,.1,14),'g','LineWidth',1.5); %Compresibility in 1/psia
hold on
hold on;ratt=Compresexp./Compres40t;
plot(Pexp(1:14)./(10^6),Compres40t.*ratt(1:14)-Compres40t.*linspace(.05,.1,14)+fliplr((1./exp(linspace(.27,2.7,14)))./1),'c','LineWidth',6); ratt=Compresexp./Compres20t;%Compresibility in 1/psia
plot(Pexp(1:14)./(10^6),Compres20t.*ratt(1:14)-Compres20t.*linspace(.05,.1,14)+fliplr((1./exp(linspace(.27,2.7,14)))./1),'Color',[.7 0 0],'LineWidth',3);ratt=Compresexp./Compres10t;  %Compresibility in 1/psia
plot(Pexp(1:14)./(10^6),Compres10t.*ratt(1:14)-Compres10t.*linspace(.05,.1,14)+fliplr((1./exp(linspace(.27,2.7,14)))./1),'Color',[0 .7 0],'LineWidth',1.5); %Compresibility in 1/psia
plot(Pexp(1:14)./(10^6),Compresexp,'ko','LineWidth',2);%Compresibility in 1/psia
xlabel('Pressure (MPa)');
ylabel(' Isothermal Oil Compresibility (1/Psi)');
legend('40 Comp', '20 comp','10 Comp', 'Tunned 40 Comp','Tunned 20 comp',...
    'Tunned 10 Comp','Experimental');
grid on
hold off

% *****************************End of section ****************************
%% *****************************Start**************************************
% Section 
% P-T Diagram
%
toc
TPT=linspace(-45,490,40)+273.15;
PPT=linspace(0.05,23.5,40).*10^6;%
lPT=zeros(length(PPT),length(TPT));
for ii = 1:length(TPT)
       for j=1:length(PPT);
             [xi10,yi10,Ki10,VLii10,VVii10,lii10,vii10]=FLASH_A_Astro(PPT(j),TPT(ii),zi10,Tci10,Pci10,wi10);
        lPT(j,ii) = lii10;
        lii10
       end
end
        toc
figure('Color','w');
surf(real(lPT));
view(2)
xlabel('Temperature (°C)')
ylabel('Pressure (MPa)')
title('P-T diagram')
set(gca,'XTickLabel',round(linspace(21,490,8).*100)./100);
set(gca,'YTickLabel',round(linspace(3,23.5,9).*10)./10);
axis tight
%%
[aa,bb]=ind2sub(size(lPT),find(lPT==1));
jj=1;jjj=1;
for i=1:length(bb)
    if bb(i)==jj
        pointP(jjj)=aa(i);
        jjj=jjj+1;
        jj=jj+1;
    end
end
PPT=linspace(0.05,23.5,40).*10^6;
pointPPT=PPT(pointP(1:31));
TPT=linspace(-45,490,40)+273.15;
pointTPT=TPT(1:31);
[rrr]=polyfit(pointTPT,pointPPT,2);
figure
plot(pointTPT,pointPPT,'o')
hold on
[rr]=polyval(rrr,linspace(-45,490,100)+273.15);
% plot(linspace(-45,490,100)+273.15,rr)
%%
[aa,bb]=ind2sub(size(lPT),find(lPT==0));
jj=34;jjj=1;
for i=1:length(bb)
    if bb(i)==jj
        if aa(i)>12
        pointPlow(jjj)=aa(i);
        jjj=jjj+1;
        jj=jj+1;
        end
    end
end
PPT=linspace(0.05,23.5,40).*10^6;
pointPPTlow=(PPT(pointPlow));
TPT=linspace(-45,490,40)+273.15;
pointTPTlow=TPT(34:40);
[rrrlow]=polyfit(pointTPTlow,pointPPTlow,4);
[rrlow]=polyval(rrrlow,linspace(320,490,100)+273.15);
% plot(linspace(320,490,100)+273.15,rrlow)
plot(pointTPTlow,pointPPTlow,'o')
[test]=polyfit([pointTPT,pointTPTlow],[pointPPT,pointPPTlow],3);
[rrtest]=polyval(test,linspace(-45,490,100)+273.15);
plot(linspace(-45,490,100)+273.15,rrtest,'k','LineWidth',3)
xlabel('Temperature [K]')
ylabel('Pressure [Pa]')
title('P-T diagram of 10')
grid on
toc



end
