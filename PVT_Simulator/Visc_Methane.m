function [mu_C1] = ViscMethane(Tref,rho_ref)

% disp('Correlation of Hanley et al. for methane')
% computes reference viscosity for CS using the correlation of Hanley et
% al. Cyrogenics, July 1975
%Tref=400.; %K
%rho_ref=0.0000406*300/400; %g/cm^3
rho_c=16.043/99.2; %g/cm^3
% parameters for dilute gas coefficient
GV=[-209097.5 264726.9 -147281.8 47167.40 -9491.872 1219.979 -96.27993 4.274152 -0.08141531];
% compute dilute gas coefficient
visc0=0.;
for i=1:9
    exp1=-1+(i-1)*1./3.;
    visc0=visc0+GV(i)*Tref^exp1;
end

% parameters for first density correction term
Avisc1=1.696985927;
Bvisc1=-0.133372346;
Cvisc1=1.4;
Fvisc1=168.0;
% first density coefficient
visc1 = Avisc1+Bvisc1*(Cvisc1-log(Tref/Fvisc1))^2.;

% parameters for the viscosity remainder
j1=-10.35060586;
j2=17.571599671;
j3=-3019.3918656;
j4=188.73011594;
j5=0.042903609488;
j6=145.29023444;
j7=6127.6818706;
% viscosity remainder
theta=(rho_ref-rho_c)/rho_c;
visc2 = rho_ref^0.1*(j2+j3/Tref^1.5)+theta*rho_ref^0.5*(j5+j6/Tref+j7/Tref^2.);
visc2= exp(visc2)-1.;
visc2=exp(j1+j4/Tref)*visc2;

% methane viscosity at T and P
% multiply by 10-4 to convert to cP
mu_C1 = (visc0+visc1+visc2)*1e-4;

end

