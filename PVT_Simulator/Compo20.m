%% *****************************Start**************************************

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
Mwi20(20) = sum(zi(20:202).*Mwi(20:202))/sum(zi(20:202));

VL20 = zeros(1,length(Pexp));
VV20 = zeros(1,length(Pexp));
l20 = zeros(1,length(Pexp));
v20 = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi20,yi20,Ki20,VLii20,VVii20,lii20,vii20]=FLASH_A_Astro(Pexp(i),T,zi20,Tci20,Pci20,wi20);
    VL20(i) = VLii20;
    VV20(i) = VVii20;
    l20(i) = lii20;
    v20(i) = 1-l20(i);
end
VL20 = VL20.*l20;
VV20 = VV20.*v20;
VT20 = VL20+VV20;
Vb20 = VT20(14);
RelV20 = VT20./Vb20;

Compres20 = (((1./VL20(1:14)).*diff(VL20(1:15)))./6892.8571).*10^6;
save('Comp20.mat')
