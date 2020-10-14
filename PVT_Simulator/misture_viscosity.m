 function [Viscmix] = Viscmix(Tci,Pci,Vmix,zi,Mwi,wi,T,P)

DenTcmix=0;
NumTcmix=0;
for ii = 1:length(zi)
    for jj = 1:length(zi)
        NumTcmix=NumTcmix+((zi(ii)*zi(jj))*((Tci(ii)*Tci(jj))^.5)*...
            (((Tci(ii)/Pci(ii))^(1/3)+(Tci(jj)/Pci(jj))^(1/3))^3));
        DenTcmix=DenTcmix+(zi(ii)*zi(jj))*((((Tci(ii)/Pci(ii))^(1/3)+...
            (Tci(jj)/Pci(jj))^(1/3))^3));
    end
end
Tcmix=NumTcmix/DenTcmix;%K
Pc_mix = 8*(NumTcmix)/(((DenTcmix)^2));%Pa

Mn_avg = sum(zi.*Mwi); 
Mw_avg = sum(zi.*(Mwi.^2))/Mn_avg;
Mw_mix = (1.304*10^-4)*(Mw_avg^2.303-Mn_avg^2.303)+Mn_avg; %g/mole

% Methane properties
Vc_ref = 99.2;%cm^3/mol
Mw_ref = 16.4; %g/mole
Pc_ref = 46;%bar
Pc_ref = Pc_ref*10^5;%Pa
Tc_ref = 190.4;%K
w_ref = 0.011;%adim
z_ref = 1; %adim

Tr = T*Tc_ref/Tcmix;%K
Pr = P*Pc_ref/Pc_mix;%Pa


[VL]=FLASH_single_component(Pr,Tr,z_ref,Tc_ref,Pc_ref,w_ref); %Input just methane properties
VLr = VL*1000000;% to convert m^3/mole to cm3/mole
rhoR = (Vc_ref/VLr); %Adim

% Volume correction should be done before alfa_mix calculation

alfa_mix = 1+(7.378*10^-3)*rhoR^1.847*Mw_mix^0.5173;%g/mole
alfa0 = 1+0.031*(rhoR^1.847);%afim

T0 = Tr*alfa0/alfa_mix;%K mole/g
rho0 = (Vmix*Mw_ref*Mw_ref/(VLr*Vc_ref*Mw_mix))*(alfa0/alfa_mix);%mole/cm^3 

% [mu_C1] = ViscMethane(Tref,rho_ref); it gives the vicoaity in cP
[mu_ref] = ViscMethane(T0,rho0);%cP
mu_ref = mu_ref *100;
Viscmix =((Tcmix/Tc_ref)^(-1/6))*((Pc_mix/Pc_ref)^(2/3))*sqrt(Mw_mix/Mw_ref)*(alfa_mix/alfa0)*mu_ref;%cP
