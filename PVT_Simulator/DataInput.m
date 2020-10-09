
%% *****************************Start**************************************
% Section
% Wellstream data

% 1    2   3   4   5   6   7   8   9   10  11  12
% N2, CO2, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10
Weight_fraction = [0.357 0.101 6.173 3.688 4.479 3.546 2.721 4.235 4.993 4.557 3.564 3.549];%percentage
Weight_fraction = Weight_fraction./100;%fraction
Mole_fraction = [1.150 0.208 34.735 11.072 9.169 5.508 3.405 4.490 4.560 3.720 2.540 2.265 17.178];%percentage
Mole_fraction = Mole_fraction./100;%fraction
Mw = [28.013 44.010 16.043 30.070 44.097 58.123 72.150 85.147 98.836 110.573 126.647 141.450]; %g/mole
Density = [0.8069 0.8172 0.3 0.3563 0.5072 0.58 0.6283 0.7109 0.7332 0.7502 0.7535 0.7684];%g/cm^3
T = 162;%F
T = (T-32)/1.8+273.15;%K

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
