%% *****************************Start**************************************
% Section 
% Tuning  plus fraction component
%
% For 20 components
wi20(20)=wi20(20)*2;
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
save('Comp20t.mat')
