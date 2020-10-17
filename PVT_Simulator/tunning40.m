%% *****************************Start**************************************
% Section 
% Tuning  plus fraction component
%
%For 40 components
wi40(40)=wi40(40).*2;
VL40t = zeros(1,length(Pexp));
VV40t = zeros(1,length(Pexp));
l40t = zeros(1,length(Pexp));
v40t = zeros(1,length(Pexp));
for i = 1:length(Pexp)
    [xi40t,yi40t,Ki40t,VLii40t,VVii40t,lii40t,vii40t]=FLASH_A_tuned(Pexp(i),T,zi40,Tci40,Pci40,wi40);
    VL40t(i) = VLii40t;
    VV40t(i) = VVii40t;
    l40t(i) = lii40t;
    v40t(i) = 1-l40t(i);
end
VL40t = VL40t.*l40t;
VV40t = VV40t.*v40t;
VT40t = VL40t+VV40t;
Vb40t = VT40t(14);
RelV40t = VT40t./Vb40t;
Compres40t = (((1./VL40t(1:14)).*diff(VL40t(1:15)))./6892.8571).*10^6;
save('Comp40t.mat')
