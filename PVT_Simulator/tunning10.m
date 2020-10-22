%% *****************************Start**************************************
% Section 
% Tuning  plus fraction component
%
%For 10 components
load('Comp10.mat');

wi10(10)=wi10(10).*2;
VL10t = zeros(1,length(Pexp));
VV10t = zeros(1,length(Pexp));
l10t = zeros(1,length(Pexp));
v10t = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi10t,yi10t,Ki10t,VLii10t,VVii10t,lii10t,vii10t]=FLASH_A_tuned(Pexp(i),T,zi10,Tci10,Pci10,wi10);
    VL10t(i) = VLii10t;
    VV10t(i) = VVii10t;
    l10t(i) = lii10t;
    v10t(i) = 1-l10t(i);
end
VL10t = VL10.*l10t;
VV10t = VV10t.*v10t;
VT10t = VL10t+VV10t;
Vb10t = VT10t(14);
RelV10t = VT10t./Vb10t;
Compres10t = (((1./VL10t(2:15)).*diff(VL10t(1:15)))./6892.8571).*10^6;
save('Comp10t.mat')