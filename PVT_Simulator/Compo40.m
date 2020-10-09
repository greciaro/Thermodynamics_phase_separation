% *****************************Start**************************************
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
Mwi40(40) = sum(zi(40:202).*Mwi(40:202))/sum(zi(40:202));
VL40 = zeros(1,length(Pexp));
VV40 = zeros(1,length(Pexp));
l40 = zeros(1,length(Pexp));
v40 = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi40,yi40,Ki40,VLii40,VVii40,lii40,vii40]=FLASH_A_Astro(Pexp(i),T,zi40,Tci40,Pci40,wi40);
    VL40(i) = VLii40;
    VV40(i) = VVii40;
    l40(i) = lii40;
    v40(i) = 1-l40(i);
end
VL40 = VL40.*l40;
VV40 = VV40.*v40;
VT40 = VL40+VV40;
Vb40 = VT40(14);
RelV40 = VT40./Vb40;


Compres40 = (((1./VL40(1:14)).*diff(VL40(1:15)))./6892.8571).*10^6;
save('Comp40.mat')
