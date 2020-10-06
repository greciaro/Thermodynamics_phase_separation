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
Mwi10(4) = sum(zi(4:5).*Mwi(4:5))/zi10(4);

Tci10(5) = sum(zi(6:8).*Tci(6:8))/sum(zi(6:8));
Pci10(5) = sum(zi(6:8).*Pci(6:8))/sum(zi(6:8));
wi10(5) = sum(zi(6:8).*wi(6:8))/sum(zi(6:8));
zi10(5) = sum(zi(6:8));
Mwi10(5) = sum(zi(6:8).*Mwi(6:8))/zi10(5);

Tci10(6) = sum(zi(9:11).*Tci(9:11))/sum(zi(9:11));
Pci10(6) = sum(zi(9:11).*Pci(9:11))/sum(zi(9:11));
wi10(6) = sum(zi(9:11).*wi(9:11))/sum(zi(9:11));
zi10(6) = sum(zi(9:11));
Mwi10(6) = sum(zi(9:11).*Mwi(9:11))/zi10(6);

Tci10(7) = sum(zi(12:15).*Mwi(12:15).*Tci(12:15))/sum(zi(12:15).*Mwi(12:15));
Pci10(7) = sum(zi(12:15).*Mwi(12:15).*Pci(12:15))/sum(zi(12:15).*Mwi(12:15));
wi10(7) = sum(zi(12:15).*Mwi(12:15).*wi(12:15))/sum(zi(12:15).*Mwi(12:15));
zi10(7) = sum(zi(12:15));
Mwi10(7) = sum(zi(12:15).*Mwi(12:15))/zi10(7);

Tci10(8) = sum(zi(16:20).*Mwi(16:20).*Tci(16:20))/sum(zi(16:20).*Mwi(16:20));
Pci10(8) = sum(zi(16:20).*Mwi(16:20).*Pci(16:20))/sum(zi(16:20).*Mwi(16:20));
wi10(8) = sum(zi(16:20).*Mwi(16:20).*wi(16:20))/sum(zi(16:20).*Mwi(16:20));
zi10(8) = sum(zi(16:20));
Mwi10(8) = sum(zi(12:15).*Mwi(12:15))/zi10(8);

Tci10(9) = sum(zi(21:28).*Mwi(21:28).*Tci(21:28))/sum(zi(21:28).*Mwi(21:28));
Pci10(9) = sum(zi(21:28).*Mwi(21:28).*Pci(21:28))/sum(zi(21:28).*Mwi(21:28));
wi10(9) = sum(zi(21:28).*Mwi(21:28).*wi(21:28))/sum(zi(21:28).*Mwi(21:28));
zi10(9) = sum(zi(21:28));
Mwi10(9) = sum(zi(21:28).*Mwi(21:28))/zi10(9);

Tci10(10) = sum(zi(29:202).*Mwi(29:202).*Tci(29:202))/sum(zi(29:202).*Mwi(29:202));
Pci10(10) = sum(zi(29:202).*Mwi(29:202).*Pci(29:202))/sum(zi(29:202).*Mwi(29:202));
wi10(10) = sum(zi(29:202).*Mwi(29:202).*wi(29:202))/sum(zi(29:202).*Mwi(29:202));
zi10(10) = sum(zi(29:202));
Mwi10(10) = sum(zi(29:202).*Mwi(29:202))/zi10(10);

T = 162;%F
T = (T-32)/1.8+273.15;%K

VL10 = zeros(1,length(Pexp));
VV10 = zeros(1,length(Pexp));
l10 = zeros(1,length(Pexp));
v10 = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi10,yi10,Ki10,VLii10,VVii10,lii10,vii10]=FLASH_A_Astro(Pexp(i),T,zi10,Tci10,Pci10,wi10);
    VL10(i) = VLii10;
    VV10(i) = VVii10;
    l10(i) = lii10;
    v10(i) = 1-l10(i);
end
VL10 = VL10.*l10;
VV10 = VV10.*v10;
VT10 = VL10+VV10;
Vb10 = VT10(14);
RelV10 = VT10./Vb10;

Compres10 = (((1./VL10(2:15)).*diff(VL10(1:15))./(6892.8571))).*10^6;
save('Comp10.mat');
